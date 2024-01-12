package worker

import (
	"fmt"
	"github.com/Taichidasheen/read_predict/pkg/opt"
	"github.com/biogo/hts/sam"
	tf "github.com/wamuir/graft/tensorflow"
	"log"
	"regexp"
	"time"
)

/*
Freeze with bug for ML tag,2024-01-04
	INPUT: Topbam: the original unaligned HiFi Bam file with fi/fp/ri/rp signals
	OUTPUT: modTopbam: modification called bam file tagged by MM,ML. And the fi/fp/ri/rp signals can be kept or removed according to KiFlag.
	Compared with def Aligned_HiFiReadsMeth, the input Topbam file is unmapped/unaligned HiFi reads.
	So, we don't need to decode CIGAR string to locate CpG positions on a read, but using python string function
*/

type TopHiFiPredictWorker struct {
	model      *tf.SavedModel
	record     *sam.Record
	resultChan chan *sam.Record
	//mappingQ   int
	//minSubDep  int
	//maxSubDep  int
	//radius     int
	//scaleFlag  bool
	//keepK      string
	opts opt.Options
	err  error
}

func NewTopHiFiPredictWorker(model *tf.SavedModel, record *sam.Record, resultChan chan *sam.Record, opts opt.Options) TopHiFiPredictWorker {
	return TopHiFiPredictWorker{
		model:      model,
		record:     record,
		resultChan: resultChan,
		opts:       opts,
		//mappingQ:   mappingQ,
		//minSubDep:  minSubDep,
		//maxSubDep:  maxSubDep,
		//radius:     radius,
		//scaleFlag:  scaleFlag,
		//keepK:      keepK,
	}
}

func (w *TopHiFiPredictWorker) Task(num int) {
	model := w.model
	record := w.record
	radius := w.opts.Radius
	scaleFlag := w.opts.ScaleFlag
	recordTag, err := extractRecordTag(record)
	if err != nil {
		log.Printf("extractRecordTag err:%v", err)
		w.err = err
		return
	}
	fn := recordTag.Fn
	rn := recordTag.Rn
	totalSubreadsDep := fn + rn
	if totalSubreadsDep >= int32(w.opts.MinSubDep) && totalSubreadsDep <= int32(w.opts.MaxSubDep) && fn >= 1 && rn >= 1 {
		readFiList := recordTag.Fi
		readFpList := recordTag.Fp
		readRiList := recordTag.Ri
		readRpList := recordTag.Rp
		readSeqList := record.Seq.Expand()
		readIsReverse := false
		readQueryLength := len(readSeqList)

		cpgPosOnRead := findCpGPos(string(readSeqList), radius)
		//log.Printf("readName:%s, len(cpgPosOnRead):%d", record.Name, len(cpgPosOnRead))
		if len(cpgPosOnRead) >= 1 {
			var xReads [][][]float32
			//start := time.Now()
			for _, posOnRead := range cpgPosOnRead {

				//remove cpg at the head or tail of a read, which causing out of range index
				if posOnRead < radius+5 {
					log.Printf("posOnRead head remove, readname:%s, posOnRead:%d", record.Name, posOnRead)
					continue
				}
				if posOnRead > readQueryLength-radius-5 {
					log.Printf("posOnRead tail remove, readname:%s, posOnRead:%d", record.Name, posOnRead)
					continue
				}
				//log.Printf("readName:%s, posOnRead:%d", readName, posOnRead)

				feature, err := HiFiRead_cpg_K_Feature(posOnRead, readIsReverse, radius, readQueryLength, readSeqList, readFiList, readFpList, readRiList, readRpList, scaleFlag)
				if err != nil {
					log.Printf("HiFiRead_cpg_K_Feature err:%v, read name:%s", err, record.Name)
					continue
				}
				//log.Printf("record name:%s, xxxx feature:%+v", record.Name, feature)
				//mat := formatClosedZMW(feature, float32(totalSubreadsDep))
				//if len(mat) == 13 {
				//	xReads = append(xReads, mat)
				//}
				mat := formatTransposedClosedZMW(feature, float32(totalSubreadsDep))
				xReads = append(xReads, mat)
			}
			//input := transpose3D(xReads)
			//transpose3D(xReads)
			//log.Printf("readName:%s, constructPredictInput cost:%v", record.Name, time.Since(start))
			log.Printf("readName:%s,  len(cpgPosOnRead):%d, len(xReads):%d", record.Name, len(cpgPosOnRead), len(xReads))
			//predict
			probes, err := predict(model, xReads)
			if err != nil {
				log.Printf("predict err:%+v, read name:%s, input:%+v", err, record.Name, xReads)
				w.err = err
				return
			}
			//计算ML tag和MM tag
			mlTag, err := genMLTag(probes)
			if err != nil {
				log.Printf("genMLTag err:%+v, read name:%s", err, record.Name)
				w.err = err
				return
			}
			mmTag, err := genMMTag(string(readSeqList))
			if err != nil {
				log.Printf("genMMTag err:%+v, read name:%s", err, record.Name)
				w.err = err
				return
			}
			record.AuxFields = append(record.AuxFields, mlTag)
			record.AuxFields = append(record.AuxFields, mmTag)
			if w.opts.KeepK == "remove" {
				removeRecordTag(record)
			}
			//写出到channel
			w.resultChan <- record
		}
	}
}

func removeRecordTag(record *sam.Record) {
	auxes := record.AuxFields
	var newAuxes sam.AuxFields
	for _, aux := range auxes {
		tag := aux.Tag()
		if tag[0] == 'f' && tag[1] == 'i' {
			continue
		}
		if tag[0] == 'f' && tag[1] == 'p' {
			continue
		}
		if tag[0] == 'r' && tag[1] == 'i' {
			continue
		}
		if tag[0] == 'r' && tag[1] == 'p' {
			continue
		}
		newAuxes = append(newAuxes, aux)
	}
	record.AuxFields = newAuxes
}

// 速度比较慢，耗时差不多100～200微秒
func findCpGPosByRegexp(readSeqString string, radius int) []int {

	start := time.Now()
	defer func() {
		log.Println("findCpGPosByRegexp cost:", time.Since(start))
	}()

	cpgString := "CG"
	readSeqLen := len(readSeqString)
	re := regexp.MustCompile(cpgString)
	matches := re.FindAllStringIndex(readSeqString, -1)

	var cpgPosOnRead []int
	for _, match := range matches {
		if match[0] < radius {
			log.Printf("fitlered match, index:%d, readSeqLen:%d", match[0], readSeqLen)
			continue
		}
		if match[0] > (readSeqLen - radius) {
			log.Printf("fitlered match, index:%d, readSeqLen:%d", match[0], readSeqLen)
			continue
		}

		cpgPosOnRead = append(cpgPosOnRead, match[0])
	}
	return cpgPosOnRead
}

// 速度快一些，耗时20微秒左右
func findCpGPos(readSeqString string, radius int) []int {

	/*start := time.Now()
	defer func() {
		log.Println("findCpGPos cost:", time.Since(start))
	}()*/

	cpgString := "CG"
	readSeqLen := len(readSeqString)

	begin := radius
	end := readSeqLen - radius

	var cpgPosOnRead []int

	for i := begin; i <= end; i++ {
		if readSeqString[i:i+2] == cpgString {
			cpgPosOnRead = append(cpgPosOnRead, i)
		}
	}
	return cpgPosOnRead
}

// 构造 13 x winsize矩阵
func formatClosedZMW(feature *Feature, npasses float32) [][]float32 {

	/*start := time.Now()
	defer func() {
		log.Println("formatClosedZMW cost:", time.Since(start))
	}()*/

	var matrix [][]float32
	fTemplateSeq := feature.TemplateSeq
	fIPDList := feature.TemplateIPDList
	fPWList := feature.TemplatePWList
	rTemplateSeq := feature.ComTemplateSeq
	rIPDList := feature.ComTemplateIPDList
	rPWList := feature.ComTemplatePWList

	winsize := len(fIPDList)
	if len(fTemplateSeq) != winsize || len(rTemplateSeq) != winsize {
		log.Printf("winsize:%d, fTemplateSeq:%s, rTemplateSeq:%s", winsize, string(fTemplateSeq), string(rTemplateSeq))
		return nil
	}
	seqDict := map[byte]int{'A': 0, 'T': 1, 'C': 2, 'G': 3}

	fseqMatArray := make([][]float32, 4)
	for i := 0; i < 4; i++ {
		fseqMatArray[i] = make([]float32, winsize)
	}

	rseqMatArray := make([][]float32, 4)
	for i := 0; i < 4; i++ {
		rseqMatArray[i] = make([]float32, winsize)
	}

	var npassesArr []float32

	//构造fseqMatArray，rseqMatArray，npassesArr
	for i := 0; i < winsize; i++ {
		fbase := fTemplateSeq[i]
		fseqMatArray[seqDict[fbase]][i] = 1
		rbase := rTemplateSeq[i]
		rseqMatArray[seqDict[rbase]][i] = 1
		npassesArr = append(npassesArr, npasses)
	}

	//13 x winsize matrix
	matrix = append(matrix, fseqMatArray...)
	matrix = append(matrix, fIPDList, fPWList)
	matrix = append(matrix, rseqMatArray...)
	matrix = append(matrix, rIPDList, rPWList, npassesArr)

	return matrix
}

func transpose2D(plane [][]float32) [][]float32 {
	rows := len(plane)
	cols := len(plane[0])

	result := make([][]float32, cols)
	for i := range result {
		result[i] = make([]float32, rows)
	}

	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			result[j][i] = plane[i][j]
		}
	}
	return result
}

// 实现 np.transpose(X_reads, (0, 2, 1))
func transpose3D(cube [][][]float32) [][][]float32 {
	/*start := time.Now()
	defer func() {
		log.Println("transpose3D cost:", time.Since(start))
	}()*/

	var result [][][]float32
	for _, plane := range cube {
		posedPlane := transpose2D(plane)
		result = append(result, posedPlane)
	}
	return result
}

// 构造 winsize x 13 矩阵
func formatTransposedClosedZMW(feature *Feature, npasses float32) [][]float32 {

	/*start := time.Now()
	defer func() {
		log.Println("formatTransposedClosedZMW cost:", time.Since(start))
	}()*/

	fTemplateSeq := feature.TemplateSeq
	fIPDList := feature.TemplateIPDList
	fPWList := feature.TemplatePWList
	rTemplateSeq := feature.ComTemplateSeq
	rIPDList := feature.ComTemplateIPDList
	rPWList := feature.ComTemplatePWList

	winsize := len(fIPDList)
	if len(fTemplateSeq) != winsize || len(rTemplateSeq) != winsize {
		log.Printf("winsize:%d, fTemplateSeq:%s, rTemplateSeq:%s", winsize, string(fTemplateSeq), string(rTemplateSeq))
		return nil
	}
	seqDict := map[byte]int{'A': 0, 'T': 1, 'C': 2, 'G': 3}

	//winsize x 13 matrix
	matrix := make([][]float32, winsize)
	for i, _ := range matrix {
		matrix[i] = make([]float32, 13)
	}

	for i := 0; i < winsize; i++ {
		//ATCG, fi, fp
		fbase := fTemplateSeq[i]
		matrix[i][seqDict[fbase]] = 1
		matrix[i][4] = fIPDList[i]
		matrix[i][5] = fPWList[i]

		//ATCG, ri, rp
		rbase := rTemplateSeq[i]
		matrix[i][6+seqDict[rbase]] = 1
		matrix[i][10] = rIPDList[i]
		matrix[i][11] = rPWList[i]
		matrix[i][12] = npasses
	}

	return matrix
}

func predict(model *tf.SavedModel, inputData interface{}) ([]float32, error) {
	//log.Printf("model:%+v", model)
	//log.Printf("inputData:%v", inputData)
	start := time.Now()
	defer func() {
		log.Println("predict cost:", time.Since(start))
	}()

	tensor, err := tf.NewTensor(inputData)
	if err != nil {
		log.Printf("tf.NewTensor err: %v", err)
		return nil, err
	}

	result, err := model.Session.Run(
		map[tf.Output]*tf.Tensor{
			// python版tensorflow/keras中定义的输入层input_layer
			model.Graph.Operation("serving_default_input_1").Output(0): tensor,
		},
		[]tf.Output{
			// python版tensorflow/keras中定义的输出层output_layer
			model.Graph.Operation("StatefulPartitionedCall").Output(0),
		},
		nil,
	)
	if err != nil {
		log.Printf("model.Session.Run err:%v", err)
		return nil, err
	}
	//log.Printf("result[0].Value:%+v", result[0].Value())
	if len(result) < 1 {
		log.Printf("predict get empty result:%+v, input:%+v", result, inputData)
		return nil, fmt.Errorf("predict get empty result")
	}
	scores, ok := result[0].Value().([][]float32)
	if !ok {
		log.Printf("not expected format, result:%+v, input:%v", result, inputData)
		return nil, fmt.Errorf("not expected format")
	}
	probes := make([]float32, len(scores))
	for i, arr := range scores {
		probes[i] = arr[0]
	}
	return probes, nil
}

func batchPredict(model *tf.SavedModel, inputs []interface{}) (interface{}, error) {
	return nil, nil
}
