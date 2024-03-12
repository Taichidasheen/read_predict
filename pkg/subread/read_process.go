package subread

import (
	"fmt"
	"github.com/Taichidasheen/read_predict/pkg/opt"
	"github.com/Taichidasheen/read_predict/pkg/pool"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
	tf "github.com/wamuir/graft/tensorflow"
	"io"
	"log"
	"os"
	"sync"
)

func ReadSubreadsBam(wg *sync.WaitGroup, recordChan chan *sam.Record, bamReader *bam.Reader, topN int) {
	defer wg.Done()
	defer close(recordChan)

	var count int
	fmt.Printf("begin count...")
	for {
		record, err := bamReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Printf("error reading bam: %v", err)
			return
		}
		count++
		if count%10000 == 0 {
			log.Printf("count:%d", count)
		}
		if topN != 0 && count == topN+1 {
			log.Printf("has processed topN:%d rows, break...", topN)
			break
		}
		recordChan <- record
	}
}

func ProcessFeatureRecordChan(wg *sync.WaitGroup, resultChan chan string, recordChan chan *sam.Record,
	opts opt.Options, refOrderMap map[string]int, cgListMap map[string][]int) {
	defer wg.Done()
	defer close(resultChan)

	concurrency := opts.Processor

	var processedCnt int

	bufferSize := 4000
	var bufferCnt int

	bufferedRecords := make(chan *sam.Record, bufferSize)
	recordPositions := make(chan *RecordPosition, bufferSize)
	recordPositionDict := &RecordPositionDict{
		ChrMap:               make(map[string]interface{}),
		CpgLocationPositions: make(map[string][]*LocatedPosition),
	}

	for record := range recordChan {
		bufferedRecords <- record
		bufferCnt++
		processedCnt++

		if bufferCnt == bufferSize {
			var wg4Cpg sync.WaitGroup
			wg4Cpg.Add(concurrency)
			for i := 0; i < concurrency; i++ {
				go func() {
					defer wg4Cpg.Done()
					for re := range bufferedRecords {
						recordPos := FindLocatedCpgs(re, opts, cgListMap)
						if recordPos != nil {
							recordPositions <- recordPos
						}
					}
				}()
			}
			close(bufferedRecords)
			wg4Cpg.Wait()
			close(recordPositions)

			log.Printf("processRecordPositions start, processed count:%d, len(recordPositions):%d", processedCnt, len(recordPositions))
			//看看是否有需要写出的部分
			topDict, err := processRecordPositions(recordPositionDict, recordPositions, refOrderMap)
			if err != nil {
				log.Printf("processRecordPositions err:%v", err)
				return
			}

			//重新初始化bufferedRecords和recordPositions
			bufferedRecords = make(chan *sam.Record, bufferSize)
			recordPositions = make(chan *RecordPosition, bufferSize)
			//重置bufferCnt
			bufferCnt = 0

			if topDict == nil {
				log.Printf("no data can output now")
				continue
			}
			var wg4Out sync.WaitGroup
			wg4Out.Add(2)
			cpgOutputChan := make(chan *CpgOutput, 2000)
			//生成zmwLine
			go processAlignedSubreadsFeature(&wg4Out, topDict, cpgOutputChan, opts)
			//sort并写到cpgOutputChan
			go sortCpgOutput2Chan(&wg4Out, cpgOutputChan, resultChan, refOrderMap)

			wg4Out.Wait()
			log.Printf("processRecordPositions end, processed count:%d", processedCnt)

		}

	}

	//已经读到文件末尾或处理到topN
	if len(bufferedRecords) > 0 {

		var wg4Cpg sync.WaitGroup
		wg4Cpg.Add(concurrency)
		for i := 0; i < concurrency; i++ {
			go func() {
				defer wg4Cpg.Done()
				for re := range bufferedRecords {
					recordPos := FindLocatedCpgs(re, opts, cgListMap)
					if recordPos != nil {
						recordPositions <- recordPos
					}
				}
			}()
		}
		close(bufferedRecords)
		wg4Cpg.Wait()
		close(recordPositions)

		//最后一次output不用重新初始化bufferedRecords和recordPositions

		log.Println("final output...")

		topDict, err := processFinalRecordPositions(recordPositionDict, recordPositions, refOrderMap)
		if err != nil {
			log.Printf("processRecordPositions err:%v", err)
			return
		}
		if topDict == nil {
			log.Printf("no data can output now")
			return
		}
		var wg4Fin sync.WaitGroup
		wg4Fin.Add(2)
		cpgOutputChan := make(chan *CpgOutput, 2000)
		//生成zmwLine
		go processAlignedSubreadsFeature(&wg4Fin, topDict, cpgOutputChan, opts)
		//sort并写到resultChan
		go sortCpgOutput2Chan(&wg4Fin, cpgOutputChan, resultChan, refOrderMap)
		wg4Fin.Wait()
	}

}

func ProcessPredictRecordChan(wg *sync.WaitGroup, resultChan chan string, recordChan chan *sam.Record,
	opts opt.Options, refOrderMap map[string]int, cgListMap map[string][]int, closedModel, openModel *tf.SavedModel) {
	defer wg.Done()
	defer close(resultChan)

	concurrency := opts.Processor

	var processedCnt int

	bufferSize := 4000
	var bufferCnt int

	bufferedRecords := make(chan *sam.Record, bufferSize)
	recordPositions := make(chan *RecordPosition, bufferSize)
	recordPositionDict := &RecordPositionDict{
		ChrMap:               make(map[string]interface{}),
		CpgLocationPositions: make(map[string][]*LocatedPosition),
	}

	for record := range recordChan {
		bufferedRecords <- record
		bufferCnt++
		processedCnt++

		if bufferCnt == bufferSize {
			var wg4Cpg sync.WaitGroup
			wg4Cpg.Add(concurrency)
			for i := 0; i < concurrency; i++ {
				go func() {
					defer wg4Cpg.Done()
					for re := range bufferedRecords {
						recordPos := FindLocatedCpgs(re, opts, cgListMap)
						if recordPos != nil {
							recordPositions <- recordPos
						}
					}
				}()
			}
			close(bufferedRecords)
			wg4Cpg.Wait()
			close(recordPositions)

			log.Printf("processRecordPositions start, processed count:%d, len(recordPositions):%d", processedCnt, len(recordPositions))
			//看看是否有需要写出的部分
			topDict, err := processRecordPositions(recordPositionDict, recordPositions, refOrderMap)
			if err != nil {
				log.Printf("processRecordPositions err:%v", err)
				return
			}

			//重新初始化bufferedRecords和recordPositions
			bufferedRecords = make(chan *sam.Record, bufferSize)
			recordPositions = make(chan *RecordPosition, bufferSize)
			//重置bufferCnt
			bufferCnt = 0

			if topDict == nil {
				log.Printf("no data can output now")
				continue
			}
			var wg4Out sync.WaitGroup
			wg4Out.Add(2)
			cpgOutputChan := make(chan *CpgOutput, 2000)
			//生成zmwLine
			go processAlignedSubreadsPredict(&wg4Out, topDict, cpgOutputChan, opts, closedModel, openModel)
			//sort并写到cpgOutputChan
			go sortCpgOutput2Chan(&wg4Out, cpgOutputChan, resultChan, refOrderMap)

			wg4Out.Wait()
			log.Printf("processRecordPositions end, processed count:%d", processedCnt)

		}

	}

	//已经读到文件末尾或处理到topN
	if len(bufferedRecords) > 0 {

		var wg4Cpg sync.WaitGroup
		wg4Cpg.Add(concurrency)
		for i := 0; i < concurrency; i++ {
			go func() {
				defer wg4Cpg.Done()
				for re := range bufferedRecords {
					recordPos := FindLocatedCpgs(re, opts, cgListMap)
					if recordPos != nil {
						recordPositions <- recordPos
					}
				}
			}()
		}
		close(bufferedRecords)
		wg4Cpg.Wait()
		close(recordPositions)

		//最后一次output不用重新初始化bufferedRecords和recordPositions

		log.Println("final output...")

		topDict, err := processFinalRecordPositions(recordPositionDict, recordPositions, refOrderMap)
		if err != nil {
			log.Printf("processRecordPositions err:%v", err)
			return
		}
		if topDict == nil {
			log.Printf("no data can output now")
			return
		}
		var wg4Fin sync.WaitGroup
		wg4Fin.Add(2)
		cpgOutputChan := make(chan *CpgOutput, 2000)
		//生成zmwLine
		go processAlignedSubreadsPredict(&wg4Fin, topDict, cpgOutputChan, opts, closedModel, openModel)
		//sort并写到resultChan
		go sortCpgOutput2Chan(&wg4Fin, cpgOutputChan, resultChan, refOrderMap)
		wg4Fin.Wait()
	}

}

func ReadSubreadsBamAndProcess(wg *sync.WaitGroup, resultChan chan string,
	opts opt.Options, bamFilePath string, topN int, cgListMap map[string][]int) {
	defer wg.Done()
	defer close(resultChan)

	//bamFilePath := "/storage/yangjianLab/caoyujie/project/meth/Bam_file/4_hifi_subreads/HG01109_WT_hifi/5x.sort.HG01109_WT_hifi.aln.bam"
	bamFile, err := os.Open(bamFilePath)
	if err != nil {
		log.Printf("could not open file %q:", err)
		return
	}
	defer bamFile.Close()
	ok, err := bgzf.HasEOF(bamFile)
	if err != nil {
		log.Printf("could not open file %q:", err)
		return
	}
	if !ok {
		log.Printf("file has no bgzf magic block: may be truncated")
		return
	}
	//读写并发度
	//concurrency := 8
	concurrency := 5
	//bam reader
	bamReader, err := bam.NewReader(bamFile, concurrency)
	if err != nil {
		log.Printf("could not read bam:%v", err)
		return
	}
	defer bamReader.Close()

	//获取ref
	bamHeader := bamReader.Header()
	refs := bamHeader.Refs()
	refOrderMap := make(map[string]int)
	for i := 0; i < len(refs); i++ {
		refChr := refs[i].Name()
		refOrderMap[refChr] = i
	}
	//log.Println("refOrderMap:", refOrderMap)

	var count int
	fmt.Printf("begin count...")

	bufferSize := 4000

	bufferedRecords := make(chan *sam.Record, bufferSize)
	recordPositions := make(chan *RecordPosition, bufferSize)
	recordPositionDict := &RecordPositionDict{
		ChrMap:               make(map[string]interface{}),
		CpgLocationPositions: make(map[string][]*LocatedPosition),
	}

	for {
		record, err := bamReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Printf("error reading bam: %v", err)
			return
		}
		count++
		if count%1000 == 0 {
			log.Printf("count:%d", count)
		}
		if topN != 0 && count == topN+1 {
			log.Printf("has processed topN:%d rows, break...", topN)
			break
		}
		bufferedRecords <- record
		//bufferSize个record处理一次
		if len(bufferedRecords) == bufferSize {
			var wg4Cpg sync.WaitGroup
			wg4Cpg.Add(10)
			for i := 0; i < 10; i++ {
				go func() {
					defer wg4Cpg.Done()
					for re := range bufferedRecords {
						recordPos := FindLocatedCpgs(re, opts, cgListMap)
						if recordPos != nil {
							recordPositions <- recordPos
						}
					}
				}()
			}
			close(bufferedRecords)
			wg4Cpg.Wait()
			close(recordPositions)

			log.Printf("processRecordPositions start, count:%d, len(recordPositions):%d", count, len(recordPositions))
			//看看是否有需要写出的部分
			topDict, err := processRecordPositions(recordPositionDict, recordPositions, refOrderMap)
			if err != nil {
				log.Printf("processRecordPositions err:%v", err)
				return
			}

			//重新初始化bufferedRecords和recordPositions
			bufferedRecords = make(chan *sam.Record, bufferSize)
			recordPositions = make(chan *RecordPosition, bufferSize)

			if topDict == nil {
				log.Printf("no data can output now")
				continue
			}
			var wg4Out sync.WaitGroup
			wg4Out.Add(2)
			cpgOutputChan := make(chan *CpgOutput, 2000)
			//生成zmwLine
			go processAlignedSubreadsFeature(&wg4Out, topDict, cpgOutputChan, opts)
			//sort并写到cpgOutputChan
			go sortCpgOutput2Chan(&wg4Out, cpgOutputChan, resultChan, refOrderMap)

			wg4Out.Wait()
			log.Printf("processRecordPositions end, count:%d", count)

		}

	}
	//已经读到文件末尾或处理到topN
	if len(bufferedRecords) > 0 {

		var wg4Cpg sync.WaitGroup
		wg4Cpg.Add(10)
		for i := 0; i < 10; i++ {
			go func() {
				defer wg4Cpg.Done()
				for re := range bufferedRecords {
					recordPos := FindLocatedCpgs(re, opts, cgListMap)
					if recordPos != nil {
						recordPositions <- recordPos
					}
				}
			}()
		}
		close(bufferedRecords)
		wg4Cpg.Wait()
		close(recordPositions)

		//最后一次output不用重新初始化bufferedRecords和recordPositions

		log.Println("final output...")

		topDict, err := processFinalRecordPositions(recordPositionDict, recordPositions, refOrderMap)
		if err != nil {
			log.Printf("processRecordPositions err:%v", err)
			return
		}
		if topDict == nil {
			log.Printf("no data can output now")
			return
		}
		var wg4Fin sync.WaitGroup
		wg4Fin.Add(2)
		cpgOutputChan := make(chan *CpgOutput, 2000)
		//生成zmwLine
		go processAlignedSubreadsFeature(&wg4Fin, topDict, cpgOutputChan, opts)
		//sort并写到resultChan
		go sortCpgOutput2Chan(&wg4Fin, cpgOutputChan, resultChan, refOrderMap)
		wg4Fin.Wait()
	}

}

func processAlignedSubreadsFeature(wg *sync.WaitGroup, topDict *TopRecordPositionDict,
	cpgOutputChan chan *CpgOutput, opts opt.Options) {
	defer wg.Done()
	//记得close resultChan, 否则会deadlock
	defer close(cpgOutputChan)

	concurrency := opts.Processor

	pool := pool.New(concurrency)

	for _, chrCpg := range topDict.ChrCpgs {
		locatedPositions := topDict.CpgLocationPositions[chrCpg]
		w := NewAlignedSubreadsFeatureWorker(locatedPositions, cpgOutputChan, opts)
		pool.Run(&w)
	}
	pool.Shutdown()

}

func processAlignedSubreadsPredict(wg *sync.WaitGroup, topDict *TopRecordPositionDict,
	cpgOutputChan chan *CpgOutput, opts opt.Options, closedModel, openModel *tf.SavedModel) {
	defer wg.Done()
	//记得close resultChan, 否则会deadlock
	defer close(cpgOutputChan)

	concurrency := opts.Processor

	pool := pool.New(concurrency)

	for _, chrCpg := range topDict.ChrCpgs {
		locatedPositions := topDict.CpgLocationPositions[chrCpg]
		w := NewAlignedSubreadsPredictWorker(closedModel, openModel, locatedPositions, cpgOutputChan, opts)
		pool.Run(&w)
	}
	pool.Shutdown()

}
