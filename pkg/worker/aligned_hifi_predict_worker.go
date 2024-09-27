package worker

import (
	"fmt"
	"github.com/Taichidasheen/read_predict/pkg/cpgpos"
	"github.com/Taichidasheen/read_predict/pkg/feature"
	"github.com/Taichidasheen/read_predict/pkg/opt"
	"github.com/Taichidasheen/read_predict/pkg/predict"
	"github.com/Taichidasheen/read_predict/pkg/record_flag"
	"github.com/Taichidasheen/read_predict/pkg/record_tag"
	"github.com/Taichidasheen/read_predict/pkg/util"
	"github.com/biogo/hts/sam"
	"github.com/rs/zerolog/log"
	tf "github.com/wamuir/graft/tensorflow"
	"strings"
)

type AlignedHiFiPredictWorker struct {
	model          *tf.SavedModel
	record         *sam.Record
	cgListMap      map[string][]int
	textResultChan chan string
	bamResultChan  chan *sam.Record
	opts           opt.Options
	err            error
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
	}
}

func (w *AlignedHiFiPredictWorker) Task(num int) {
	model := w.model
	record := w.record
	radius := w.opts.Radius
	scaleFlag := w.opts.ScaleFlag
	keepK := w.opts.KeepK

	predictFlag := false //记录是否发生了predict动作

	if !record_flag.IsSecondary(record.Flags) && !record_flag.IsSupplementary(record.Flags) && int(record.MapQ) > w.opts.MappingQ && record_flag.MatchingRatio(record) >= 0.85 {
		recordTag, err := record_tag.ExtractRecordTag(record)
		if err != nil {
			log.Error().Msgf("extractRecordTag err:%v", err)
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
			log.Debug().Msgf("alnRefStart:%d, alnRefEnd:%d", alnRefStart, alnRefEnd)
			cgList := w.cgListMap[alnRefChr]
			overlappingCpg := cpgpos.FindOverlappingCpg(cgList, alnRefStart, alnRefEnd)
			log.Debug().Msgf("overlappingCpg:%+v", overlappingCpg)
			if len(overlappingCpg) >= 1 {
				readFiList := recordTag.Fi
				readFpList := recordTag.Fp
				readRiList := recordTag.Ri
				readRpList := recordTag.Rp
				readSeqList := record.Seq.Expand()
				readIsReverse := record_flag.IsReverse(record.Flags)
				readQueryLength := len(readSeqList)
				readCigar := record.Cigar
				readName := record.Name

				locatedCpgs, cpgPosOnSeq := cpgpos.LocateCpgPosOnSeq(alnRefStart, readCigar, overlappingCpg)
				log.Debug().Msgf("readName:%s, len(overlappingCpg):%d, len(locatedCpgs):%d", readName, len(overlappingCpg), len(locatedCpgs))
				log.Debug().Msgf("locatedCpgs:%+v", locatedCpgs)

				if len(locatedCpgs) >= 1 {
					var xReads [][][]float32
					var featurePosOnSeq []int //pos on the seq
					var featurePosOnRef []int //pos on the reference
					for _, cpg := range locatedCpgs {
						//heading or tailing removing
						posOnSeq := cpgPosOnSeq[cpg]
						if posOnSeq < radius+5 {
							log.Warn().Msgf("posOnRead heading removing, readname:%s, posOnSeq:%d", readName, posOnSeq)
							continue
						}
						if posOnSeq > readQueryLength-radius-5 {
							log.Warn().Msgf("posOnRead heading removing, readname:%s, posOnSeq:%d", readName, posOnSeq)
							continue
						}
						//log.Printf("readName:%s, cpg:%d, posOnRead:%d", readName, cpg, posOnRead)

						feat, err := feature.HiFiRead_cpg_K_Feature(posOnSeq, readIsReverse, radius, readQueryLength, readSeqList, readFiList, readFpList, readRiList, readRpList, scaleFlag)
						if err != nil {
							log.Error().Msgf("HiFiRead_cpg_K_Feature err:%v, read name:%s", err, record.Name)
							continue
						}

						mat := feature.FormatTransposedClosedZMW(feat, float32(totalSubreadsDep))
						xReads = append(xReads, mat)
						featurePosOnSeq = append(featurePosOnSeq, posOnSeq)
						featurePosOnRef = append(featurePosOnRef, cpg)
					}

					log.Debug().Msgf("readName:%s,  len(cpgPosOnSeq):%d, len(xReads):%d", record.Name, len(cpgPosOnSeq), len(xReads))

					predictFlag = true //修改标记为需要进行predict
					//predict
					probes, err := predict.Predict(model, xReads)
					if err != nil {
						log.Error().Msgf("predict err:%+v, read name:%s, input:%+v", err, record.Name, xReads)
						w.err = err
						return
					}

					//output
					if w.opts.OutputType == "MoleculeLevel" {
						moleculeLevelOut(record, recordTag, featurePosOnRef, probes, w.textResultChan)
					}
					if w.opts.OutputType == "ModBam" {
						err = modBamOut(record, keepK, featurePosOnSeq, probes, w.bamResultChan)
						if err != nil {
							log.Error().Msgf("modBamOut err:%+v, read name:%s", err, record.Name)
							w.err = err
							return
						}
					}
				}
			}
		}
	}
	//如果当前record没有进行predict，则去掉ML，MM tag后原样输出
	if predictFlag == false && w.opts.OutputType == "ModBam" {
		log.Debug().Msgf("no predict process,direct output, record name:%s", record.Name)
		//如果原来的文件里有MM和ML，先去掉
		record_tag.RemoveMMMLTag(record)
		if w.opts.KeepK == "remove" {
			record_tag.RemoveRecordTag(record)
		}
		//写出到channel
		w.bamResultChan <- record
	}

}

func moleculeLevelOut(record *sam.Record, recordTag *record_tag.RecordTag, featurePosOnRef []int, probes []float32, textResultChan chan string) {
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

	allLines := make([]string, len(featurePosOnRef))

	for i, pos := range featurePosOnRef {
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
	readIsReverse := record_flag.IsReverse(record.Flags)
	//readQueryLength := len(readSeqList)

	var topStrand []byte
	//var featurePosOnRead []int

	// 2024-01-05: need to distinguish +- strand, since MM,ML tag is with respected to top_strand(5'->3')
	if readIsReverse {
		probes = util.ReverseSlice(probes)
		topStrand = seqComplementary(util.ReverseSliceByte(readSeqList))
		//featurePosOnRead = convertPosOnSeq2PosOnReadSlice(readQueryLength, featurePosOnSeq)
	} else {
		topStrand = readSeqList
		//featurePosOnRead = featurePosOnSeq
	}

	//计算ML tag和MM tag
	mlTag, err := record_tag.GenMLTag(probes)
	if err != nil {
		log.Error().Msgf("genMLTag err:%+v, read name:%s", err, record.Name)
		return err
	}
	mmTag, err := record_tag.GenMMTag(string(topStrand))
	if err != nil {
		log.Error().Msgf("genMMTag err:%+v, read name:%s", err, record.Name)
		return err
	}
	//如果原来的文件里有MM和ML，先去掉
	record_tag.RemoveMMMLTag(record)

	record.AuxFields = append(record.AuxFields, mlTag)
	record.AuxFields = append(record.AuxFields, mmTag)
	if keepK == "remove" {
		record_tag.RemoveRecordTag(record)
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

// 计算反向互补序列
func seqComplementary(seq []byte) []byte {
	complementedSeq := make([]byte, len(seq))
	for i := 0; i < len(seq); i++ {
		switch seq[i] {
		case 'A':
			complementedSeq[i] = 'T'
		case 'T':
			complementedSeq[i] = 'A'
		case 'C':
			complementedSeq[i] = 'G'
		case 'G':
			complementedSeq[i] = 'C'
		default:
			log.Warn().Msgf("invalid bytes:%v", seq[i])
		}
	}
	return complementedSeq
}
