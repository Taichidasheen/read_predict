package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"github.com/Taichidasheen/read_predict/pkg/opt"
	"github.com/Taichidasheen/read_predict/pkg/pool"
	"github.com/Taichidasheen/read_predict/pkg/subread"
	"github.com/Taichidasheen/read_predict/pkg/worker"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
	tf "github.com/wamuir/graft/tensorflow"
	"io"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

func main() {

	var inType string
	var bamFilePath string
	var outPrefix string
	var cpgBed string
	var minSubreadsDepth int
	var maxSubreadsDepth int
	var radius int
	var mappingQuality int
	var baseQuality int
	var chromosome string
	var processor int
	var scale bool
	var datatype string
	var modelJson string
	var modelWeights string
	var modelDir string

	var outputType string
	var keepK string

	//debug param
	var cpuprofile string
	var topN int

	//input related parameters
	flag.StringVar(&inType, "inType", "", "aligned or unaligned")
	flag.StringVar(&bamFilePath, "bam", "", "aligned/unaligned BAM file with kinetic signals, fi/fp/ri/rp for hifi reads, and ip/pw for subreads")
	flag.StringVar(&cpgBed, "cg", "", "CpG position, column1 is chr, column2 is pos, such as chr1 132. Here pos is 1-based coordinate on Rreference's forward strand")
	flag.StringVar(&datatype, "type", "SUBREAD", "data mode, SUBREAD or HIFI")
	flag.IntVar(&radius, "r", 10, "cpg +/- r")
	flag.IntVar(&minSubreadsDepth, "min", 1, "min subreads depth (by-strand) for a ZMW to be included")
	flag.IntVar(&maxSubreadsDepth, "max", 60, "max subreads depth (by-strand) for a ZMW to be excluded")
	flag.IntVar(&mappingQuality, "maq", 30, "min subreads mapping quality from aligner")
	flag.IntVar(&baseQuality, "bpq", 0, "min base quality required")
	flag.StringVar(&chromosome, "chr", "whole_genome", "processing chromosome")
	flag.StringVar(&modelJson, "j", "", "json file of the model, json file")
	flag.StringVar(&modelWeights, "w", "", "weights of the model, h5 file")
	flag.StringVar(&modelDir, "model", "", "model dir")

	// processing parameters
	flag.IntVar(&processor, "p", 0, "Parallelism processors")

	// output related parameters
	flag.StringVar(&outPrefix, "o", "", "[*outprefix*].modification.bam for ModBam outputType.  "+
		"[*outprefix*].Kmat.txt.gz for Feature outputType. [*outprefix*].SingleMol.pre.txt.gz for MoleculeLevel outputType")
	flag.BoolVar(&scale, "sc", false, "flag to indicate scale the signal of each subreads, (x-mean)/std. Without this tag is raw signal value")
	flag.StringVar(&outputType, "oe", "", "Output Types: ModBam,Feature,MoleculeLevel. ModBam: Modification Bam file. | Feature: Feature matrix. | MoleculeLevel: Molecule level modification prediction")
	flag.StringVar(&keepK, "keepK", "remove", "remove or keep. flag to indicate keep or remove the kinetic signals in the output Bam file")

	//debug parameters
	flag.StringVar(&cpuprofile, "cpuprofile", "", "write cpu profile to this file")
	flag.IntVar(&topN, "topN", 0, "just process top N rows")

	flag.Parse()

	fmt.Println("inType:", inType)
	fmt.Println("bamFilePath:", bamFilePath)
	fmt.Println("cpgBed:", cpgBed)
	fmt.Println("minSubreadsDepth:", minSubreadsDepth)
	fmt.Println("maxSubreadsDepth:", maxSubreadsDepth)
	fmt.Println("radius:", radius)
	fmt.Println("mappingQuality:", mappingQuality)
	fmt.Println("baseQuality:", baseQuality)
	fmt.Println("chromosome:", chromosome)
	fmt.Println("processor:", processor)
	fmt.Println("datatype:", datatype)
	fmt.Println("modelDir:", modelDir)

	fmt.Println("outPrefix:", outPrefix)
	fmt.Println("scale:", scale)
	fmt.Println("outputType:", outputType)
	fmt.Println("keepK:", keepK)

	fmt.Println("cpuprofile:", cpuprofile)
	fmt.Println("topN:", topN)

	fmt.Println("maxProcs:", runtime.GOMAXPROCS(0))

	if cpuprofile != "" {
		f, err := os.Create(cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	opts := opt.Options{
		MappingQ:   mappingQuality,
		MinSubDep:  minSubreadsDepth,
		MaxSubDep:  maxSubreadsDepth,
		Radius:     radius,
		ScaleFlag:  scale,
		KeepK:      keepK,
		OutputType: outputType,
		Processor:  processor,
	}

	var wg sync.WaitGroup

	var taskName string
	if datatype == "HIFI" {
		if inType == "unaligned" {
			if outputType == "Feature" {
				taskName = inType + "_" + outputType
				wg.Add(3)
				recordChan := make(chan *sam.Record, 1000)
				resultChan := make(chan string, 3000)
				//predictResultChan := make(chan *sam.Record, 1000)

				//读取bam文件
				go readBam(&wg, recordChan, bamFilePath, topN)

				//异步处理record
				go processTopHiFiFeature(&wg, recordChan, resultChan, opts)

				//resultPath := "/storage/yangjianLab/westlakechat/subreads_locate/result.txt"
				//异步写入result
				go writeTextResult(&wg, outPrefix, resultChan)
				//bamHeader := bamReader.Header()
				//go writeBamRecord(&wg, bamHeader, outPrefix, predictResultChan)
			}
		}

		if inType == "aligned" {
			if outputType == "Feature" {
				taskName = inType + "_" + outputType
				//load cgList
				cgListMap, err := buildCgListMap(cpgBed, chromosome)
				if err != nil {
					log.Fatalf("buildCgListMap err:%v", err)
					return
				}
				log.Println("len(cgListMap):", len(cgListMap))

				wg.Add(3)
				recordChan := make(chan *sam.Record, 1000)
				resultChan := make(chan string, 3000)
				//predictResultChan := make(chan *sam.Record, 1000)

				//读取bam文件
				go readBam(&wg, recordChan, bamFilePath, topN)

				//异步处理record
				go processAlignedHiFiFeature(&wg, cgListMap, recordChan, resultChan, opts)

				//resultPath := "/storage/yangjianLab/westlakechat/subreads_locate/result.txt"
				//异步写入result
				go writeTextResult(&wg, outPrefix, resultChan)
				//bamHeader := bamReader.Header()
				//go writeBamRecord(&wg, bamHeader, outPrefix, predictResultChan)
			}
		}
	}

	if datatype == "SUBREAD" {
		if outputType == "Feature" {
			taskName = datatype + "_" + inType + "_" + outputType
			//load cgList
			cgListMap, err := buildCgListMap(cpgBed, chromosome)
			if err != nil {
				log.Fatalf("buildCgListMap err:%v", err)
				return
			}
			log.Println("len(cgListMap):", len(cgListMap))

			bamFile, err := os.Open(bamFilePath)
			if err != nil {
				log.Fatalf("could not open file %q:", err)
				return
			}
			defer bamFile.Close()
			ok, err := bgzf.HasEOF(bamFile)
			if err != nil {
				log.Fatalf("could not open file %q:", err)
				return
			}
			if !ok {
				log.Printf("file has no bgzf magic block: may be truncated")
				return
			}
			//读写并发度
			concurrency := 5
			//bam reader
			bamReader, err := bam.NewReader(bamFile, concurrency)
			if err != nil {
				log.Fatalf("could not read bam:%v", err)
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

			wg.Add(3)
			recordChan := make(chan *sam.Record, 1000)
			resultChan := make(chan string, 3000)

			//读取bam文件
			go subread.ReadSubreadsBam(&wg, recordChan, bamReader, topN)
			//go subread.ReadSubreadsBamAndProcess(&wg, resultChan, opts, bamFilePath, topN, cgListMap)
			//处理record
			go subread.ProcessRecordChan(&wg, resultChan, recordChan, opts, refOrderMap, cgListMap)

			//异步写入result
			go writeTextResult(&wg, outPrefix, resultChan)

		}

		//predict
		if outputType == "MoleculeLevel" {

		}
	}

	//load模型
	/*model, err := loadModel(modelDir, []string{"serve"})
	if err != nil {
		log.Printf("loadModel err:%v, modelDir:%s", err, modelDir)
	}
	log.Printf("loaded model:%v", model)*/

	if taskName == "" {
		log.Println("unsupported params, no task processed")
		return
	} else {
		log.Printf("begin to process task:%s", taskName)
	}

	wg.Wait()

	fmt.Println("exit...")

}

func loadModel(modelPath string, modelNames []string) (*tf.SavedModel, error) {
	model, err := tf.LoadSavedModel(modelPath, modelNames, nil) // 载入模型
	if err != nil {
		log.Printf("LoadSavedModel err: %v", err)
		return nil, err
	}

	log.Println("list possible ops in graphs")
	for _, op := range model.Graph.Operations() {
		log.Printf("Op name: %v", op.Name())
	}

	return model, nil
}

func buildCgListMap(cgFile string, processingChr string) (map[string][]int, error) {
	start := time.Now()
	defer func() {
		log.Println("buildCgListMap cost:", time.Since(start))
	}()

	cgListMap := make(map[string][]int)
	// 打开文件
	file, err := os.Open(cgFile)
	if err != nil {
		log.Println("Error opening file:", err)
		return nil, err
	}
	defer file.Close()

	// 创建一个 Scanner 来读取文件内容
	scanner := bufio.NewScanner(file)

	// 逐行读取文件内容
	for scanner.Scan() {
		line := scanner.Text()
		parts := strings.Split(line, "\t")
		if len(parts) != 2 {
			log.Printf("invalid line:%s", line)
			continue
		}
		pos, err := strconv.Atoi(parts[1])
		if err != nil {
			log.Printf("invalid line:%s", line)
			continue
		}
		chr := parts[0]
		if processingChr == "whole_genome" || chr == processingChr {
			cgListMap[chr] = append(cgListMap[chr], pos)
		}
	}

	if err = scanner.Err(); err != nil {
		log.Println("Error reading file:", err)
		return nil, err
	}
	//sort
	for chr, _ := range cgListMap {
		sort.Ints(cgListMap[chr])
	}
	return cgListMap, nil
}

func readBam(wg *sync.WaitGroup, recordChan chan *sam.Record, bamFilePath string, topN int) {
	defer wg.Done()
	defer close(recordChan)

	//bamFilePath := "/storage/yangjianLab/caoyujie/project/meth/Bam_file/4_hifi_subreads/HG01109_WT_hifi/5x.sort.HG01109_WT_hifi.aln.bam"
	bamFile, err := os.Open(bamFilePath)
	if err != nil {
		log.Fatalf("could not open file %q:", err)
		return
	}
	defer bamFile.Close()
	ok, err := bgzf.HasEOF(bamFile)
	if err != nil {
		log.Fatalf("could not open file %q:", err)
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
		log.Fatalf("could not read bam:%v", err)
		return
	}
	defer bamReader.Close()

	var count int
	fmt.Printf("begin count...")
	for {
		record, err := bamReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("error reading bam: %v", err)
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

func writeTextResult(wg *sync.WaitGroup, outPrefix string, resultChan chan string) {
	defer wg.Done()

	file, err := os.Create(outPrefix)
	if err != nil {
		log.Fatalf("could not open file %q:", err)
		return
	}
	defer file.Close()

	// 创建一个写入器
	writer := bufio.NewWriter(file)
	for line := range resultChan {
		_, err := writer.WriteString(line + "\n")
		if err != nil {
			log.Println("Error writing line:", err)
			return
		}
	}
	// 将缓冲区的数据刷新到文件中
	err = writer.Flush()
	if err != nil {
		log.Println("Error flushing writer:", err)
		return
	}
}

func writeGzipResult(wg *sync.WaitGroup, outPrefix string, resultChan chan string) {
	defer wg.Done()

	resultPath := outPrefix + ".Kmat.txt.gz"

	file, err := os.Create(resultPath)
	if err != nil {
		log.Fatalf("could not open file %q:", err)
		return
	}
	defer file.Close()

	// 创建一个写入器
	gzipWriter := gzip.NewWriter(file)
	defer gzipWriter.Close()

	for line := range resultChan {
		_, err := gzipWriter.Write([]byte(line + "\n"))
		if err != nil {
			log.Println("Error writing line:", err)
			return
		}
	}
	// 将缓冲区的数据刷新到文件中
	err = gzipWriter.Flush()
	if err != nil {
		log.Println("Error flushing writer:", err)
		return
	}
}

func writeBamRecord(wg *sync.WaitGroup, bamHeader *sam.Header, outBamPath string, predictResultChan chan *sam.Record) {
	defer wg.Done()
	outBam, err := os.OpenFile(outBamPath, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0644)
	if err != nil {
		log.Fatalf("could not open file %q:", err)
		return
	}
	defer outBam.Close()

	//header, _ := sam.NewHeader(nil, nil)
	bamWriter, err := bam.NewWriter(outBam, bamHeader, 5)
	if err != nil {
		log.Fatalf("could not write bam:%v", err)
		return
	}
	defer bamWriter.Close()

	for record := range predictResultChan {
		err := bamWriter.Write(record)
		if err != nil {
			log.Printf("write bam err:%+v, record:%s", err, record.Name)
			continue
		}
	}
}

func processAlignedHiFiFeature(wg *sync.WaitGroup, cgListMap map[string][]int, recordChan chan *sam.Record,
	resultChan chan string, opts opt.Options) {
	defer wg.Done()

	concurrency := opts.Processor

	pool := pool.New(concurrency)

	for record := range recordChan {
		w := worker.NewAlignedHiFiFeatureWorker(record, cgListMap, resultChan, opts)
		pool.Run(&w)
	}
	pool.Shutdown()
	//记得close resultChan, 否则会deadlock
	close(resultChan)
}

func processTopHiFiFeature(wg *sync.WaitGroup, recordChan chan *sam.Record, resultChan chan string, opts opt.Options) {
	defer wg.Done()

	concurrency := opts.Processor

	pool := pool.New(concurrency)
	for record := range recordChan {
		w := worker.NewTopHiFiFeatureWorker(record, resultChan, opts)
		pool.Run(&w)
	}
	pool.Shutdown()
	//记得close resultChan, 否则会deadlock
	close(resultChan)
}

func processTopHiFiPredict(wg *sync.WaitGroup, model *tf.SavedModel, recordChan chan *sam.Record,
	resultChan chan *sam.Record, opts opt.Options) {
	defer wg.Done()

	concurrency := opts.Processor

	pool := pool.New(concurrency)
	for record := range recordChan {
		w := worker.NewTopHiFiPredictWorker(model, record, resultChan, opts)
		pool.Run(&w)
	}
	pool.Shutdown()
	//记得close resultChan, 否则会deadlock
	close(resultChan)
}
