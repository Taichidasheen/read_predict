package main

import (
	"bufio"
	"flag"
	"fmt"
	"github.com/Taichidasheen/subreads_locate/pkg/pool"
	"github.com/Taichidasheen/subreads_locate/pkg/worker"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
	"io"
	"log"
	"os"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

func main() {

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

	flag.StringVar(&bamFilePath, "bam", "", "aligned BAM file")
	flag.StringVar(&outPrefix, "o", "", "")
	flag.StringVar(&cpgBed, "cg", "", "CpG position, chr pos, pos is 1-based coordinate on reference's forward strand")
	flag.IntVar(&minSubreadsDepth, "min", 1, "min subreads depth (by-strand) for a ZMW to be included")
	flag.IntVar(&maxSubreadsDepth, "max", 60, "max subreads depth (by-strand) for a ZMW to be excluded")
	flag.IntVar(&radius, "r", 10, "radius of the study window")
	flag.IntVar(&mappingQuality, "maq", 30, "min subreads mapping quality from aligner")
	flag.IntVar(&baseQuality, "bpq", 0, "min base quality required")
	flag.StringVar(&chromosome, "chr", "whole_genome", "processing chromosome")
	flag.IntVar(&processor, "p", 0, "Parallelism processors")
	flag.BoolVar(&scale, "sc", false, "flag to indicate scale the signal of each subreads, (x-mean)/std. Without this tag is raw signal value")
	flag.StringVar(&datatype, "type", "SUBREAD", "data mode, SUBREAD or HIFI")
	flag.Parse()

	fmt.Println("bamFilePath:", bamFilePath)
	fmt.Println("outPrefix:", outPrefix)
	fmt.Println("cpgBed:", cpgBed)
	fmt.Println("minSubreadsDepth:", minSubreadsDepth)
	fmt.Println("maxSubreadsDepth:", maxSubreadsDepth)
	fmt.Println("radius:", radius)
	fmt.Println("mappingQuality:", mappingQuality)
	fmt.Println("baseQuality:", baseQuality)
	fmt.Println("chromosome:", chromosome)
	fmt.Println("processor:", processor)
	fmt.Println("scale:", scale)
	fmt.Println("datatype:", datatype)

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

	cgListMap, err := buildCgListMap(cpgBed, chromosome)
	if err != nil {
		log.Fatalf("buildCgListMap err:%v", err)
		return
	}

	log.Println("len(cgListMap):", len(cgListMap))

	log.Printf("maxProcs:%d", runtime.GOMAXPROCS(0))
	//读写并发度
	//concurrency := 8
	concurrency := 3

	//bam reader
	bamReader, err := bam.NewReader(bamFile, concurrency)
	if err != nil {
		log.Fatalf("could not read bam:%v", err)
		return
	}
	defer bamReader.Close()

	recordChan := make(chan *sam.Record, 500)
	resultChan := make(chan string, 3000)

	var wg sync.WaitGroup
	wg.Add(2)
	//异步处理record
	go processRecord(&wg, processor, cgListMap, recordChan, resultChan, mappingQuality, minSubreadsDepth, maxSubreadsDepth, radius, scale)

	//resultPath := "/storage/yangjianLab/westlakechat/subreads_locate/result.txt"
	//异步写入result
	go writeResult(&wg, outPrefix, resultChan)

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
		if count%100 == 0 {
			log.Printf("n:%d", count)
			break
		}
		recordChan <- record
	}
	close(recordChan)

	wg.Wait()

	fmt.Println("exit...")

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

func writeResult(wg *sync.WaitGroup, resultPath string, resultChan chan string) {
	defer wg.Done()

	file, err := os.Create(resultPath)
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

func processRecord(wg *sync.WaitGroup, concurrency int, cgListMap map[string][]int, recordChan chan *sam.Record, resultChan chan string,
	mappingQ, minSubDep, maxSubDep, radius int, scaleFlag bool) {
	defer wg.Done()

	pool := pool.New(concurrency)
	for record := range recordChan {
		w := worker.NewLocateWorker(record, cgListMap, resultChan, mappingQ, minSubDep, maxSubDep, radius, scaleFlag)
		pool.Run(&w)
	}
	pool.Shutdown()
	//记得close resultChan, 否则会deadlock
	close(resultChan)
}
