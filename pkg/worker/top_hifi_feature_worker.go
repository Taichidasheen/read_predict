package worker

import (
	"fmt"
	"github.com/Taichidasheen/read_predict/pkg/feature"
	"github.com/Taichidasheen/read_predict/pkg/opt"
	"github.com/Taichidasheen/read_predict/pkg/record_tag"
	"github.com/Taichidasheen/read_predict/pkg/util"
	"github.com/biogo/hts/sam"
	"github.com/rs/zerolog/log"
	"strings"
)

type TopHiFiFeatureWorker struct {
	record     *sam.Record
	resultChan chan string
	//mappingQ   int
	//minSubDep  int
	//maxSubDep  int
	//radius     int
	//scaleFlag  bool
	opts opt.Options
	err  error
}

func NewTopHiFiFeatureWorker(record *sam.Record, resultChan chan string, opts opt.Options) TopHiFiFeatureWorker {
	return TopHiFiFeatureWorker{
		record:     record,
		resultChan: resultChan,
		opts:       opts,
		//mappingQ:   mappingQ,
		//minSubDep:  minSubDep,
		//maxSubDep:  maxSubDep,
		//radius:     radius,
		//scaleFlag:  scaleFlag,
	}
}

func (w *TopHiFiFeatureWorker) Task(num int) {
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
	if totalSubreadsDep >= int32(w.opts.MinSubDep) && totalSubreadsDep <= int32(w.opts.MaxSubDep) && fn >= 1 && rn >= 1 {
		readFiList := recordTag.Fi
		readFpList := recordTag.Fp
		readRiList := recordTag.Ri
		readRpList := recordTag.Rp
		readSeqList := record.Seq.Expand()
		readIsReverse := false
		readQueryLength := len(readSeqList)
		readName := record.Name
		parts := strings.Split(readName, "/")
		zmwname := parts[0] + "_" + parts[1]

		cpgPosOnRead := findCpGPos(string(readSeqList), radius)
		log.Debug().Msgf("readName:%s, len(cpgPosOnRead):%d", record.Name, len(cpgPosOnRead))
		if len(cpgPosOnRead) >= 1 {
			for _, posOnRead := range cpgPosOnRead {
				//remove cpg at the head or tail of a read, which causing out of range index
				if posOnRead < radius+5 {
					log.Warn().Msgf("posOnRead head remove, readname:%s, posOnRead:%d", readName, posOnRead)
					continue
				}
				if posOnRead > readQueryLength-radius-5 {
					log.Warn().Msgf("posOnRead tail remove, readname:%s, posOnRead:%d", readName, posOnRead)
					continue
				}
				log.Debug().Msgf("readName:%s, posOnRead:%d", readName, posOnRead)

				feat, err := feature.HiFiRead_cpg_K_Feature(posOnRead, readIsReverse, radius, readQueryLength, readSeqList, readFiList, readFpList, readRiList, readRpList, scaleFlag)
				if err != nil {
					log.Error().Msgf("HiFiRead_cpg_K_Feature err:%v, read name:%s", err, record.Name)
					continue
				}

				//输出
				fTemplateSeqStr := string(feat.TemplateSeq)
				fIPDPart := util.FormatSlice(feat.TemplateIPDList)
				fPWPart := util.FormatSlice(feat.TemplatePWList)
				rTemplateSeqStr := string(feat.ComTemplateSeq)
				rIPDPart := util.FormatSlice(feat.ComTemplateIPDList)
				rPWPart := util.FormatSlice(feat.ComTemplatePWList)

				outputZMWLine := fmt.Sprintf("Top\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s",
					zmwname, posOnRead, fn, rn, fTemplateSeqStr, fIPDPart, fPWPart,
					rTemplateSeqStr, rIPDPart, rPWPart)
				w.resultChan <- outputZMWLine

			}
		}

	}

}
