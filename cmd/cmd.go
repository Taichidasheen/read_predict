package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"github.com/Taichidasheen/read_predict/pkg/fasta"
	"github.com/Taichidasheen/read_predict/pkg/opt"
	"github.com/Taichidasheen/read_predict/pkg/pool"
	"github.com/Taichidasheen/read_predict/pkg/subread"
	"github.com/Taichidasheen/read_predict/pkg/worker"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
	"github.com/rs/zerolog"
	"runtime"

	"github.com/rs/zerolog/log"
	tf "github.com/wamuir/graft/tensorflow"
	"io"
	"os"
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
	var chromosome string
	var processor int
	var scale bool
	var datatype string
	var hifiModelDir string
	var closedModelDir string
	var openModelDir string

	var outputType string
	var keepK string

	//debug param
	var cpuprofile string
	var topN int
	var logLevel int

	//input parameters
	flag.StringVar(&bamFilePath, "bamfile", "", "aligned/unaligned BAM file with kinetic signals, fi/fp/ri/rp for hifi reads, and ip/pw for subreads")
	flag.StringVar(&inType, "InBamStatus", "", "aligned or unaligned")
	flag.StringVar(&datatype, "readformat", "HIFI", "read type, SUBREAD or HIFI")
	flag.StringVar(&outputType, "outputType", "", "Output Types: ModBam,Feature,MoleculeLevel. ModBam: Modification Bam file. | Feature: Feature matrix. | MoleculeLevel: Molecule level modification prediction")
	flag.StringVar(&outPrefix, "outprefix", "", "[*outprefix*].modification.bam for ModBam outputType.  "+
		"[*outprefix*].Kmat.txt.gz for Feature outputType. [*outprefix*].SingleMol.pre.txt.gz for MoleculeLevel outputType")
	flag.StringVar(&keepK, "KeepKinetics", "remove", "remove or keep. flag to indicate keep or remove the kinetic signals in the output Bam file")

	//ZMWmeth model files
	flag.StringVar(&hifiModelDir, "HDir", "", "ZMWmeth-HiFi model directory")
	flag.StringVar(&closedModelDir, "CDir", "", "ZMWmeth-ClosedZMW model directory")
	flag.StringVar(&openModelDir, "ODir", "", "ZMWmeth-OpenZMW model directory")

	//control parameters
	flag.IntVar(&radius, "radius", 10, "cpg +/- r base pairs will be included in the read-level prediction")
	flag.IntVar(&minSubreadsDepth, "minsubreadsdepth", 1, "min subreads depth for a ZMW to be included")
	flag.IntVar(&maxSubreadsDepth, "maxsubreadsdepth", 60, "max subreads depth for a ZMW to be excluded")
	flag.IntVar(&mappingQuality, "MappingQuality", 30, "min mapping quality for a read to be included")
	flag.StringVar(&cpgBed, "Reference", "", "The reference genome in fasta format when the input is an aligned bam")
	flag.StringVar(&chromosome, "Chromosome", "whole_genome", "processing chromosome")
	flag.BoolVar(&scale, "Scale", false, "flag to indicate scale the signal of each subreads, (x-mean)/std. Without this tag is raw signal value")

	// processing parameters
	flag.IntVar(&processor, "Processor", 0, "Parallelism processors")

	//debug parameters
	flag.StringVar(&cpuprofile, "cpuprofile", "", "write cpu profile to this file")
	flag.IntVar(&topN, "topN", 0, "just process top N rows")
	flag.IntVar(&logLevel, "loglevel", 1, "0-debug,1-info,2-warn,3-error")

	flag.Parse()
	log.Output(os.Stdout)
	zerolog.SetGlobalLevel(zerolog.Level(logLevel))

	/*//input parameters
	fmt.Println("bamfile:", bamFilePath)
	fmt.Println("InBamStatus:", inType)
	fmt.Println("readformat:", datatype)
	fmt.Println("outputType:", outputType)
	fmt.Println("outprefix:", outPrefix)
	fmt.Println("KeepKinetics:", keepK)

	//ZMWmeth model files
	fmt.Println("HDir:", hifiModelDir)
	fmt.Println("CDir:", closedModelDir)
	fmt.Println("ODir:", openModelDir)

	//control parameters
	fmt.Println("radius:", radius)
	fmt.Println("minsubreadsdepth:", minSubreadsDepth)
	fmt.Println("maxsubreadsdepth:", maxSubreadsDepth)
	fmt.Println("MappingQuality:", mappingQuality)
	fmt.Println("Processor:", processor)
	fmt.Println("Reference:", cpgBed)
	fmt.Println("Scale:", scale)
	fmt.Println("Chromosome:", chromosome)

	fmt.Println("cpuprofile:", cpuprofile)
	fmt.Println("topN:", topN)

	fmt.Println("maxProcs:", runtime.GOMAXPROCS(0))*/

	//input parameters
	log.Info().Msgf("bamfile:%v", bamFilePath)
	log.Info().Msgf("InBamStatus:%v", inType)
	log.Info().Msgf("readformat:%v", datatype)
	log.Info().Msgf("outputType:%v", outputType)
	log.Info().Msgf("outprefix:%v", outPrefix)
	log.Info().Msgf("KeepKinetics:%v", keepK)

	//ZMWmeth model files
	log.Info().Msgf("HDir:%v", hifiModelDir)
	log.Info().Msgf("CDir:%v", closedModelDir)
	log.Info().Msgf("ODir:%v", openModelDir)

	//control parameters
	log.Info().Msgf("radius:%v", radius)
	log.Info().Msgf("minsubreadsdepth:%v", minSubreadsDepth)
	log.Info().Msgf("maxsubreadsdepth:%v", maxSubreadsDepth)
	log.Info().Msgf("MappingQuality:%v", mappingQuality)
	log.Info().Msgf("Processor:%v", processor)
	log.Info().Msgf("Reference:%v", cpgBed)
	log.Info().Msgf("Scale:%v", scale)
	log.Info().Msgf("Chromosome:%v", chromosome)

	log.Info().Msgf("cpuprofile:%v", cpuprofile)
	log.Info().Msgf("topN:%v", topN)
	log.Info().Msgf("loglevel:%v", logLevel)

	log.Info().Msgf("maxProcs:%v", runtime.GOMAXPROCS(0))

	if cpuprofile != "" {
		f, err := os.Create(cpuprofile)
		if err != nil {
			log.Err(err)
			return
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
				taskName = datatype + "_" + inType + "_" + outputType

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
				outfile := outPrefix + ".Kmat.txt"
				go writeTextResult(&wg, outfile, resultChan)
				//bamHeader := bamReader.Header()
				//go writeBamRecord(&wg, bamHeader, outPrefix, predictResultChan)
			}
			if outputType == "MoleculeLevel" {
				taskName = datatype + "_" + inType + "_" + outputType

				model, err := loadModel(hifiModelDir, []string{"serve"})
				if err != nil {
					log.Printf("loadModel err:%v, modelDir:%s", err, hifiModelDir)
					return
				}
				log.Printf("loaded model:%v", model)

				bamFile, err := os.Open(bamFilePath)
				if err != nil {
					log.Error().Msgf("could not open file %q:", err)
					return
				}
				defer bamFile.Close()
				ok, err := bgzf.HasEOF(bamFile)
				if err != nil {
					log.Error().Msgf("could not open file %q:", err)
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
					log.Error().Msgf("could not read bam:%v", err)
					return
				}
				defer bamReader.Close()

				wg.Add(3)
				recordChan := make(chan *sam.Record, 1000)
				//resultChan := make(chan string, 3000)
				predictResultChan := make(chan *sam.Record, 1000)

				//读取bam文件
				//go readBam(&wg, recordChan, bamFilePath, topN)
				go readHiFiBam(&wg, recordChan, bamReader, topN)
				//异步处理record
				go processTopHiFiPredict(&wg, model, recordChan, predictResultChan, opts)

				//resultPath := "/storage/yangjianLab/westlakechat/subreads_locate/result.txt"
				//异步写入result
				//go writeTextResult(&wg, outPrefix, resultChan)
				outFile := outPrefix + ".5mc.mod.unaligned.bam"
				bamHeader := bamReader.Header()
				go writeBamRecord(&wg, bamHeader, outFile, predictResultChan)
			}
		}

		if inType == "aligned" {
			if outputType == "Feature" {
				taskName = datatype + "_" + inType + "_" + outputType
				//load cgList
				//cgListMap, err := buildCgListMap(cpgBed, chromosome)
				cgListMap, err := buildCgListMapFromFasta(cpgBed, chromosome)
				if err != nil {
					log.Error().Msgf("buildCgListMap err:%v", err)
					return
				}
				log.Debug().Msgf("len(cgListMap):%d", len(cgListMap))

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
				outFile := outPrefix + ".Kmat.txt"
				go writeTextResult(&wg, outFile, resultChan)
				//bamHeader := bamReader.Header()
				//go writeBamRecord(&wg, bamHeader, outPrefix, predictResultChan)
			}
			if outputType == "MoleculeLevel" || outputType == "ModBam" {
				taskName = datatype + "_" + inType + "_" + outputType

				//load cgList
				//cgListMap, err := buildCgListMap(cpgBed, chromosome)
				cgListMap, err := buildCgListMapFromFasta(cpgBed, chromosome)
				if err != nil {
					log.Error().Msgf("buildCgListMap err:%v", err)
					return
				}
				log.Debug().Msgf("len(cgListMap):%d", len(cgListMap))

				model, err := loadModel(hifiModelDir, []string{"serve"})
				if err != nil {
					log.Printf("loadModel err:%v, modelDir:%s", err, hifiModelDir)
				}
				log.Printf("loaded model:%v", model)

				bamFile, err := os.Open(bamFilePath)
				if err != nil {
					log.Error().Msgf("could not open file %q:", err)
					return
				}
				defer bamFile.Close()
				ok, err := bgzf.HasEOF(bamFile)
				if err != nil {
					log.Error().Msgf("could not open file %q:", err)
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
					log.Error().Msgf("could not read bam:%v", err)
					return
				}
				defer bamReader.Close()

				wg.Add(3)
				recordChan := make(chan *sam.Record, 1000)

				predictTextResultChan := make(chan string, 3000)
				predictBamResultChan := make(chan *sam.Record, 1000)

				//读取bam文件
				//go readBam(&wg, recordChan, bamFilePath, topN)
				go readHiFiBam(&wg, recordChan, bamReader, topN)
				//异步处理record
				go processAlignedHiFiPredict(&wg, cgListMap, model, recordChan, predictTextResultChan, predictBamResultChan, opts)

				//resultPath := "/storage/yangjianLab/westlakechat/subreads_locate/result.txt"
				//异步写入result
				if outputType == "MoleculeLevel" {
					outFile := outPrefix + ".SingleMol.pre.txt"
					go writeTextResult(&wg, outFile, predictTextResultChan)
				}
				if outputType == "ModBam" {
					outFile := outPrefix + ".modification.bam"
					bamHeader := bamReader.Header()
					go writeBamRecord(&wg, bamHeader, outFile, predictBamResultChan)
				}
			}
		}
	}

	if datatype == "SUBREAD" {
		if outputType == "Feature" {
			taskName = datatype + "_" + inType + "_" + outputType
			//load cgList
			//cgListMap, err := buildCgListMap(cpgBed, chromosome)
			cgListMap, err := buildCgListMapFromFasta(cpgBed, chromosome)
			if err != nil {
				log.Error().Msgf("buildCgListMap err:%v", err)
				return
			}
			log.Debug().Msgf("len(cgListMap):%d", len(cgListMap))

			bamFile, err := os.Open(bamFilePath)
			if err != nil {
				log.Error().Msgf("could not open file %q:", err)
				return
			}
			defer bamFile.Close()
			ok, err := bgzf.HasEOF(bamFile)
			if err != nil {
				log.Error().Msgf("could not open file %q:", err)
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
				log.Error().Msgf("could not read bam:%v", err)
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
			go subread.ProcessFeatureRecordChan(&wg, resultChan, recordChan, opts, refOrderMap, cgListMap)

			//异步写入result
			outFile := outPrefix + ".Kmat.txt"
			go writeTextResult(&wg, outFile, resultChan)

		}

		//predict
		if outputType == "MoleculeLevel" {
			taskName = datatype + "_" + inType + "_" + outputType
			//load cgList
			//cgListMap, err := buildCgListMap(cpgBed, chromosome)
			cgListMap, err := buildCgListMapFromFasta(cpgBed, chromosome)
			if err != nil {
				log.Error().Msgf("buildCgListMap err:%v", err)
				return
			}
			log.Debug().Msgf("len(cgListMap):%d", len(cgListMap))

			//load model
			closedModel, err := loadModel(closedModelDir, []string{"serve"})
			if err != nil {
				log.Printf("loadModel err:%v, closedModelDir:%s", err, closedModelDir)
			}
			log.Printf("loaded closed model:%v", closedModel)

			openModel, err := loadModel(openModelDir, []string{"serve"})
			if err != nil {
				log.Printf("loadModel err:%v, openModelDir:%s", err, openModelDir)
			}
			log.Printf("loaded open model:%v", openModel)

			bamFile, err := os.Open(bamFilePath)
			if err != nil {
				log.Error().Msgf("could not open file %q:", err)
				return
			}
			defer bamFile.Close()
			ok, err := bgzf.HasEOF(bamFile)
			if err != nil {
				log.Error().Msgf("could not open file %q:", err)
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
				log.Error().Msgf("could not read bam:%v", err)
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
			go subread.ProcessPredictRecordChan(&wg, resultChan, recordChan, opts, refOrderMap, cgListMap, closedModel, openModel)

			//异步写入result
			outFile := outPrefix + ".SingleMol.pre.txt"
			go writeTextResult(&wg, outFile, resultChan)

		}
	}

	//load模型
	/*model, err := loadModel(modelDir, []string{"serve"})
	if err != nil {
		log.Printf("loadModel err:%v, modelDir:%s", err, modelDir)
	}
	log.Printf("loaded model:%v", model)*/

	if taskName == "" {
		log.Error().Msgf("unsupported params, no task processed")
		return
	} else {
		log.Printf("begin to process task:%s", taskName)
	}

	wg.Wait()

	log.Info().Msgf("exit...")

}

func loadModel(modelPath string, modelNames []string) (*tf.SavedModel, error) {
	model, err := tf.LoadSavedModel(modelPath, modelNames, nil) // 载入模型
	if err != nil {
		log.Error().Msgf("LoadSavedModel err: %v", err)
		return nil, err
	}

	log.Debug().Msgf("list possible ops in graphs")
	for _, op := range model.Graph.Operations() {
		log.Debug().Msgf("Op name: %v", op.Name())
	}

	return model, nil
}

func buildCgListMap(cgFile string, processingChr string) (map[string][]int, error) {
	start := time.Now()
	defer func() {
		log.Debug().Msgf("buildCgListMap cost:", time.Since(start))
	}()

	cgListMap := make(map[string][]int)
	// 打开文件
	file, err := os.Open(cgFile)
	if err != nil {
		log.Error().Msgf("Error opening file:%v", err)
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
		log.Error().Msgf("Error reading file:%v", err)
		return nil, err
	}
	//sort
	for chr, _ := range cgListMap {
		sort.Ints(cgListMap[chr])
	}
	return cgListMap, nil
}

func buildCgListMapFromFasta(cgFile string, processingChr string) (map[string][]int, error) {
	start := time.Now()
	defer func() {
		log.Debug().Msgf("buildCgListMapFromFasta cost:%v", time.Since(start))
	}()

	cgListMap := make(map[string][]int)

	records, err := fasta.ReadFASTA(cgFile)
	if err != nil {
		log.Printf("fasta.ReadFASTA err:%+v", err)
		return nil, err
	}

	for _, record := range records {
		chr := record.Chromosome
		if processingChr == "whole_genome" || chr == processingChr {
			// 查找 "CG" 出现的位置
			cgPositions := record.FindCGPositions()
			//fmt.Println("CG Positions:", cgPositions[:50])
			cgListMap[chr] = cgPositions
		}
	}
	return cgListMap, nil
}

func readBam(wg *sync.WaitGroup, recordChan chan *sam.Record, bamFilePath string, topN int) {
	defer wg.Done()
	defer close(recordChan)

	//bamFilePath := "/storage/yangjianLab/caoyujie/project/meth/Bam_file/4_hifi_subreads/HG01109_WT_hifi/5x.sort.HG01109_WT_hifi.aln.bam"
	bamFile, err := os.Open(bamFilePath)
	if err != nil {
		log.Error().Msgf("could not open file %q:", err)
		return
	}
	defer bamFile.Close()
	ok, err := bgzf.HasEOF(bamFile)
	if err != nil {
		log.Error().Msgf("could not open file %q:", err)
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
		log.Error().Msgf("could not read bam:%v", err)
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
			log.Error().Msgf("error reading bam: %v", err)
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

func readHiFiBam(wg *sync.WaitGroup, recordChan chan *sam.Record, bamReader *bam.Reader, topN int) {
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

func writeTextResult(wg *sync.WaitGroup, outFile string, resultChan chan string) {
	defer wg.Done()

	file, err := os.Create(outFile)
	if err != nil {
		log.Error().Msgf("could not open file %q:", err)
		return
	}
	defer file.Close()

	// 创建一个写入器
	writer := bufio.NewWriter(file)
	for line := range resultChan {
		_, err := writer.WriteString(line + "\n")
		if err != nil {
			log.Error().Msgf("Error writing line:%v", err)
			return
		}
	}
	// 将缓冲区的数据刷新到文件中
	err = writer.Flush()
	if err != nil {
		log.Error().Msgf("Error flushing writer:%v", err)
		return
	}
}

func writeGzipResult(wg *sync.WaitGroup, outPrefix string, resultChan chan string) {
	defer wg.Done()

	resultPath := outPrefix + ".Kmat.txt.gz"

	file, err := os.Create(resultPath)
	if err != nil {
		log.Error().Msgf("could not open file %q:", err)
		return
	}
	defer file.Close()

	// 创建一个写入器
	gzipWriter := gzip.NewWriter(file)
	defer gzipWriter.Close()

	for line := range resultChan {
		_, err := gzipWriter.Write([]byte(line + "\n"))
		if err != nil {
			log.Error().Msgf("Error writing line:%v", err)
			return
		}
	}
	// 将缓冲区的数据刷新到文件中
	err = gzipWriter.Flush()
	if err != nil {
		log.Error().Msgf("Error flushing writer:%v", err)
		return
	}
}

func writeBamRecord(wg *sync.WaitGroup, bamHeader *sam.Header, outBamPath string, predictResultChan chan *sam.Record) {
	defer wg.Done()
	outBam, err := os.OpenFile(outBamPath, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0644)
	if err != nil {
		log.Error().Msgf("could not open file %q:", err)
		return
	}
	defer outBam.Close()

	//header, _ := sam.NewHeader(nil, nil)
	bamWriter, err := bam.NewWriter(outBam, bamHeader, 5)
	if err != nil {
		log.Error().Msgf("could not write bam:%v", err)
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

func processAlignedHiFiPredict(wg *sync.WaitGroup, cgListMap map[string][]int, model *tf.SavedModel, recordChan chan *sam.Record,
	textResultChan chan string, bamResultChan chan *sam.Record, opts opt.Options) {
	defer wg.Done()

	concurrency := opts.Processor

	pool := pool.New(concurrency)
	for record := range recordChan {
		w := worker.NewAlignedHiFiPredictWorker(model, record, cgListMap, textResultChan, bamResultChan, opts)
		pool.Run(&w)
	}
	pool.Shutdown()
	//记得close resultChan, 否则会deadlock
	close(textResultChan)
	close(bamResultChan)
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
