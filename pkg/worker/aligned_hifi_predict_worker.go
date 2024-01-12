package worker

import (
	"fmt"
	"github.com/Taichidasheen/read_predict/pkg/opt"
	"github.com/biogo/hts/sam"
	tf "github.com/wamuir/graft/tensorflow"
	"log"
	"strings"
)

type AlignedHiFiPredictWorker struct {
	model          *tf.SavedModel
	record         *sam.Record
	cgListMap      map[string][]int
	textResultChan chan string
	bamResultChan  chan *sam.Record
	opts           opt.Options
	//mappingQ   int
	//minSubDep  int
	//maxSubDep  int
	//radius     int
	//scaleFlag  bool
	//outputType string
	err error
}

func NewAlignedHiFiPredictWorker(model *tf.SavedModel, record *sam.Record, cgListMap map[string][]int,
	textResultChan chan string, bamResultChan chan *sam.Record, opts opt.Options) AlignedHiFiPredictWorker {
	return AlignedHiFiPredictWorker{
		model:          model,
		record:         record,
		cgListMap:      cgListMap,
		textResultChan: textResultChan,
		bamResultChan:  bamResultChan,
		opts:           opts,
		//mappingQ:   mappingQ,
		//minSubDep:  minSubDep,
		//maxSubDep:  maxSubDep,
		//radius:     radius,
		//scaleFlag:  scaleFlag,
		//outputType: outputType,
	}
}

func (w *AlignedHiFiPredictWorker) Task(num int) {
	model := w.model
	record := w.record
	radius := w.opts.Radius
	scaleFlag := w.opts.ScaleFlag
	keepK := w.opts.KeepK
	if !isSecondary(record.Flags) && !isSupplementary(record.Flags) && int(record.MapQ) > w.opts.MappingQ && matchingRatio(record) >= 0.85 {
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
			alnRefChr := record.Ref.Name()
			alnRefStart := record.Pos
			alnRefEnd := record.End()
			log.Printf("alnRefStart:%d, alnRefEnd:%d", alnRefStart, alnRefEnd)
			cgList := w.cgListMap[alnRefChr]
			overlappingCpg := findOverlappingCpg(cgList, alnRefStart, alnRefEnd)
			log.Printf("overlappingCpg:%+v", overlappingCpg)
			if len(overlappingCpg) >= 1 {
				readFiList := recordTag.Fi
				readFpList := recordTag.Fp
				readRiList := recordTag.Ri
				readRpList := recordTag.Rp
				readSeqList := record.Seq.Expand()
				readIsReverse := isReverse(record.Flags)
				readQueryLength := len(readSeqList)
				readCigar := record.Cigar
				readName := record.Name

				locatedCpgs, cpgPosOnSeq := locateCpgPosOnSeq(alnRefStart, readCigar, overlappingCpg)
				log.Printf("readName:%s, len(overlappingCpg):%d, len(locatedCpgs):%d", readName, len(overlappingCpg), len(locatedCpgs))
				log.Printf("locatedCpgs:%+v", locatedCpgs)

				if len(locatedCpgs) >= 1 {
					var xReads [][][]float32
					var featurePosOnSeq []int
					for _, cpg := range locatedCpgs {
						//heading or tailing removing
						posOnSeq := cpgPosOnSeq[cpg]
						if posOnSeq < radius+5 {
							log.Printf("posOnRead heading removing, readname:%s, posOnSeq:%d", readName, posOnSeq)
							continue
						}
						if posOnSeq > readQueryLength-radius-5 {
							log.Printf("posOnRead heading removing, readname:%s, posOnSeq:%d", readName, posOnSeq)
							continue
						}
						//log.Printf("readName:%s, cpg:%d, posOnRead:%d", readName, cpg, posOnRead)

						feature, err := HiFiRead_cpg_K_Feature(posOnSeq, readIsReverse, radius, readQueryLength, readSeqList, readFiList, readFpList, readRiList, readRpList, scaleFlag)
						if err != nil {
							log.Printf("HiFiRead_cpg_K_Feature err:%v, read name:%s", err, record.Name)
							continue
						}

						mat := formatTransposedClosedZMW(feature, float32(totalSubreadsDep))
						xReads = append(xReads, mat)
						featurePosOnSeq = append(featurePosOnSeq, posOnSeq)
					}

					log.Printf("readName:%s,  len(cpgPosOnSeq):%d, len(xReads):%d", record.Name, len(cpgPosOnSeq), len(xReads))

					//predict
					probes, err := predict(model, xReads)
					if err != nil {
						log.Printf("predict err:%+v, read name:%s, input:%+v", err, record.Name, xReads)
						w.err = err
						return
					}

					//output
					if w.opts.OutputType == "MoleculeLevel" {
						moleculeLevelOut(record, recordTag, featurePosOnSeq, probes, w.textResultChan)
					}
					if w.opts.OutputType == "ModBam" {
						err = modBamOut(record, keepK, featurePosOnSeq, probes, w.bamResultChan)
						if err != nil {
							log.Printf("modBamOut err:%+v, read name:%s", err, record.Name)
							w.err = err
							return
						}
					}
				}
			}
		}
	}
}

func moleculeLevelOut(record *sam.Record, recordTag *RecordTag, featurePosOnSeq []int, probes []float32, textResultChan chan string) {
	alnRefChr := record.Ref.Name()
	readName := record.Name
	parts := strings.Split(readName, "/")
	zmwName := parts[0] + "_" + parts[1]

	fn := recordTag.Fn
	rn := recordTag.Rn
	totalSubreadsDep := fn + rn
	haplotype := recordTag.HP
	haploTypeBlock := recordTag.PS

	zmwType := "H"
	strandSign := "+"

	allLines := make([]string, len(featurePosOnSeq))

	for i, pos := range featurePosOnSeq {
		var bi int
		if probes[i] >= 0.5 {
			bi = 1
		} else {
			bi = 0
		}
		line := fmt.Sprintf("%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%v\t%d",
			alnRefChr, pos, zmwName, zmwType, totalSubreadsDep, haplotype, haploTypeBlock, strandSign, probes[i], bi)
		allLines[i] = line
	}
	allLinesStr := strings.Join(allLines, "\n")
	textResultChan <- allLinesStr
}

func modBamOut(record *sam.Record, keepK string, featurePosOnSeq []int, probes []float32, bamResultChan chan *sam.Record) error {
	readSeqList := record.Seq.Expand()
	readIsReverse := isReverse(record.Flags)
	readQueryLength := len(readSeqList)

	var topStrand []byte
	var featurePosOnRead []int

	// 2024-01-05: need to distinguish +- strand, since MM,ML tag is with respected to top_strand(5'->3')
	if readIsReverse {
		probes = reverseSlice(probes)
		topStrand = seqComplementary(reverseSliceByte(readSeqList))
		featurePosOnRead = convertPosOnSeq2PosOnReadSlice(readQueryLength, featurePosOnSeq)
	} else {
		topStrand = readSeqList
		featurePosOnRead = featurePosOnSeq
	}

	//计算ML tag和MM tag
	mlTag, err := genMLTag(probes)
	if err != nil {
		log.Printf("genMLTag err:%+v, read name:%s", err, record.Name)
		return err
	}
	mmTag, err := genAlignedMMTag(topStrand, featurePosOnRead)
	if err != nil {
		log.Printf("genMMTag err:%+v, read name:%s", err, record.Name)
		return err
	}
	record.AuxFields = append(record.AuxFields, mlTag)
	record.AuxFields = append(record.AuxFields, mmTag)
	if keepK == "remove" {
		removeRecordTag(record)
	}
	//写出到channel
	bamResultChan <- record
	return nil
}

// ConvertPosSeq2PosReadMap
// Freeze 2024-01-05
// Convert pos_on_SEQ(list) to pos_on_top_strand(5'->3')
// SEQ is from the alignment bam file
func convertPosSeq2PosReadMap(readIsReverse bool, readLength int, cpgPosOnSeq map[int]int) map[int]int {
	cpgPosOnRead := make(map[int]int)
	if readIsReverse == false {
		cpgPosOnRead = cpgPosOnSeq
	} else {
		for cpg, seqPos := range cpgPosOnSeq {
			readPos := readLength - seqPos - 1
			cpgPosOnRead[cpg] = readPos
		}
	}
	return cpgPosOnSeq
}

func convertPosOnSeq2PosOnReadSlice(readLength int, posOnSeq []int) []int {

	arrLen := len(posOnSeq)
	posOnRead := make([]int, arrLen)
	for i := arrLen - 1; i >= 0; i-- {
		posOnRead[arrLen-i-1] = readLength - posOnSeq[i] - 1
	}
	return posOnRead
}

func genAlignedMMTag(topStrand []byte, featurePosOnRead []int) (sam.Aux, error) {
	var mmArr []string

	//注意：这里是len(featurePos)-1,如果是len(featurePos)可能会数组越界
	for i := 0; i < len(featurePosOnRead)-1; i++ {
		start := featurePosOnRead[i]
		end := featurePosOnRead[i+1]
		countC := 0

		for j := start; j < end; j++ {
			if topStrand[j] == 'C' {
				countC++
			}
		}
		mmArr = append(mmArr, fmt.Sprintf("%d", countC-1))
	}
	mmVal := "C+m," + strings.Join(mmArr, ",") + ";"
	return sam.NewAux(sam.NewTag("MM"), mmVal)
}
