package worker

import (
	"github.com/Taichidasheen/read_predict/pkg/feature"
	"github.com/Taichidasheen/read_predict/pkg/opt"
	"github.com/Taichidasheen/read_predict/pkg/predict"
	"github.com/Taichidasheen/read_predict/pkg/record_tag"
	"github.com/biogo/hts/sam"
	"github.com/rs/zerolog/log"
	tf "github.com/wamuir/graft/tensorflow"

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
	opts       opt.Options
	err        error
}

func NewTopHiFiPredictWorker(model *tf.SavedModel, record *sam.Record, resultChan chan *sam.Record, opts opt.Options) TopHiFiPredictWorker {
	return TopHiFiPredictWorker{
		model:      model,
		record:     record,
		resultChan: resultChan,
		opts:       opts,
	}
}

func (w *TopHiFiPredictWorker) Task(num int) {
	model := w.model
	record := w.record
	radius := w.opts.Radius
	scaleFlag := w.opts.ScaleFlag
	recordTag, err := record_tag.ExtractRecordTag(record)
	if err != nil {
		log.Error().Msgf("extractRecordTag err:%v", err)
		w.err = err
		return
	}
	fn := recordTag.Fn
	rn := recordTag.Rn
	totalSubreadsDep := fn + rn

	predictFlag := false //记录是否发生了predict动作

	if totalSubreadsDep >= int32(w.opts.MinSubDep) && totalSubreadsDep <= int32(w.opts.MaxSubDep) && fn >= 1 && rn >= 1 {
		readFiList := recordTag.Fi
		readFpList := recordTag.Fp
		readRiList := recordTag.Ri
		readRpList := recordTag.Rp
		readSeqList := record.Seq.Expand()
		readIsReverse := false
		readQueryLength := len(readSeqList)

		cpgPosOnRead := findCpGPos(string(readSeqList), radius)
		//log.Debug().Msgf("readName:%s, len(cpgPosOnRead):%d", record.Name, len(cpgPosOnRead))
		if len(cpgPosOnRead) >= 1 {
			var xReads [][][]float32
			//start := time.Now()
			for _, posOnRead := range cpgPosOnRead {

				//remove cpg at the head or tail of a read, which causing out of range index
				if posOnRead < radius+5 {
					log.Warn().Msgf("posOnRead head remove, readname:%s, posOnRead:%d", record.Name, posOnRead)
					continue
				}
				if posOnRead > readQueryLength-radius-5 {
					log.Warn().Msgf("posOnRead tail remove, readname:%s, posOnRead:%d", record.Name, posOnRead)
					continue
				}
				//log.Debug().Msgf("readName:%s, posOnRead:%d", readName, posOnRead)

				feat, err := feature.HiFiRead_cpg_K_Feature(posOnRead, readIsReverse, radius, readQueryLength, readSeqList, readFiList, readFpList, readRiList, readRpList, scaleFlag)
				if err != nil {
					log.Error().Msgf("HiFiRead_cpg_K_Feature err:%v, read name:%s", err, record.Name)
					continue
				}
				//log.Debug().Msgf("record name:%s, xxxx feature:%+v", record.Name, feature)
				//mat := formatClosedZMW(feature, float32(totalSubreadsDep))
				//if len(mat) == 13 {
				//	xReads = append(xReads, mat)
				//}
				mat := feature.FormatTransposedClosedZMW(feat, float32(totalSubreadsDep))
				xReads = append(xReads, mat)
			}
			//input := transpose3D(xReads)
			//transpose3D(xReads)
			//log.Debug().Msgf("readName:%s, constructPredictInput cost:%v", record.Name, time.Since(start))
			log.Debug().Msgf("readName:%s,  len(cpgPosOnRead):%d, len(xReads):%d", record.Name, len(cpgPosOnRead), len(xReads))
			predictFlag = true //修改标记为需要进行predict
			//predict
			probes, err := predict.Predict(model, xReads)
			if err != nil {
				log.Error().Msgf("predict err:%+v, read name:%s, input:%+v", err, record.Name, xReads)
				w.err = err
				return
			}
			//计算ML tag和MM tag
			mlTag, err := record_tag.GenMLTag(probes)
			if err != nil {
				log.Error().Msgf("genMLTag err:%+v, read name:%s", err, record.Name)
				w.err = err
				return
			}
			mmTag, err := record_tag.GenMMTag(string(readSeqList))
			if err != nil {
				log.Error().Msgf("genMMTag err:%+v, read name:%s", err, record.Name)
				w.err = err
				return
			}
			//如果原来的文件里有MM和ML，先去掉
			record_tag.RemoveMMMLTag(record)

			record.AuxFields = append(record.AuxFields, mlTag)
			record.AuxFields = append(record.AuxFields, mmTag)
			if w.opts.KeepK == "remove" {
				record_tag.RemoveRecordTag(record)
			}
			//写出到channel
			w.resultChan <- record
		}
	}
	//如果当前record没有进行predict，则去掉ML，MM tag后原样输出
	if predictFlag == false {
		log.Debug().Msgf("no predict process,direct output, record name:%s", record.Name)
		//如果原来的文件里有MM和ML，先去掉
		record_tag.RemoveMMMLTag(record)
		if w.opts.KeepK == "remove" {
			record_tag.RemoveRecordTag(record)
		}
		//写出到channel
		w.resultChan <- record
	}
}

// 速度比较慢，耗时差不多100～200微秒
func findCpGPosByRegexp(readSeqString string, radius int) []int {

	start := time.Now()
	defer func() {
		log.Debug().Msgf("findCpGPosByRegexp cost:%v", time.Since(start))
	}()

	cpgString := "CG"
	readSeqLen := len(readSeqString)
	re := regexp.MustCompile(cpgString)
	matches := re.FindAllStringIndex(readSeqString, -1)

	var cpgPosOnRead []int
	for _, match := range matches {
		if match[0] < radius {
			log.Warn().Msgf("fitlered match, index:%d, readSeqLen:%d", match[0], readSeqLen)
			continue
		}
		if match[0] > (readSeqLen - radius) {
			log.Warn().Msgf("fitlered match, index:%d, readSeqLen:%d", match[0], readSeqLen)
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
