package worker

import (
	"fmt"
	"github.com/Taichidasheen/read_predict/pkg/cpgpos"
	"github.com/Taichidasheen/read_predict/pkg/feature"
	"github.com/Taichidasheen/read_predict/pkg/opt"
	"github.com/Taichidasheen/read_predict/pkg/record_flag"
	"github.com/Taichidasheen/read_predict/pkg/record_tag"
	"github.com/Taichidasheen/read_predict/pkg/util"
	"github.com/biogo/hts/sam"
	"log"
	"strings"
)

type AlignedHiFiFeatureWorker struct {
	record     *sam.Record
	cgListMap  map[string][]int
	resultChan chan string
	opts       opt.Options
	err        error
}

func NewAlignedHiFiFeatureWorker(record *sam.Record, cgListMap map[string][]int,
	resultChan chan string, opts opt.Options) AlignedHiFiFeatureWorker {
	return AlignedHiFiFeatureWorker{
		record:     record,
		cgListMap:  cgListMap,
		resultChan: resultChan,
		opts:       opts,
	}
}

func (w *AlignedHiFiFeatureWorker) Task(num int) {
	record := w.record
	radius := w.opts.Radius
	scaleFlag := w.opts.ScaleFlag
	if !record_flag.IsSecondary(record.Flags) && !record_flag.IsSupplementary(record.Flags) && int(record.MapQ) > w.opts.MappingQ && record_flag.MatchingRatio(record) >= 0.85 {
		recordTag, err := record_tag.ExtractRecordTag(record)
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
			overlappingCpg := cpgpos.FindOverlappingCpg(cgList, alnRefStart, alnRefEnd)
			log.Printf("overlappingCpg:%+v", overlappingCpg)
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
				parts := strings.Split(readName, "/")
				zmwname := parts[0] + "_" + parts[1]

				locatedCpgs, cpgPosOnRead := cpgpos.LocateCpgPosOnSeq(alnRefStart, readCigar, overlappingCpg)
				log.Printf("readName:%s, len(overlappingCpg):%d, len(locatedCpgs):%d", readName, len(overlappingCpg), len(locatedCpgs))
				log.Printf("locatedCpgs:%+v", locatedCpgs)

				if len(locatedCpgs) >= 1 {
					for _, cpg := range locatedCpgs {
						//heading or tailing removing
						posOnRead := cpgPosOnRead[cpg]
						if posOnRead < radius+5 {
							log.Printf("posOnRead heading removing, readname:%s, posOnRead:%d", readName, posOnRead)
							continue
						}
						if posOnRead > readQueryLength-radius-5 {
							log.Printf("posOnRead heading removing, readname:%s, posOnRead:%d", readName, posOnRead)
							continue
						}
						//log.Printf("readName:%s, cpg:%d, posOnRead:%d", readName, cpg, posOnRead)

						feat, err := feature.HiFiRead_cpg_K_Feature(posOnRead, readIsReverse, radius, readQueryLength, readSeqList, readFiList, readFpList, readRiList, readRpList, scaleFlag)
						if err != nil {
							log.Printf("HiFiRead_cpg_K_Feature err:%v, read name:%s", err, record.Name)
							continue
						}

						//输出
						fTemplateSeqStr := string(feat.TemplateSeq)
						fIPDPart := util.FormatSlice(feat.TemplateIPDList)
						fPWPart := util.FormatSlice(feat.TemplatePWList)
						rTemplateSeqStr := string(feat.ComTemplateSeq)
						rIPDPart := util.FormatSlice(feat.ComTemplateIPDList)
						rPWPart := util.FormatSlice(feat.ComTemplatePWList)

						outputZMWLine := fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s",
							alnRefChr, cpg, zmwname, fn, rn, fTemplateSeqStr, fIPDPart, fPWPart,
							rTemplateSeqStr, rIPDPart, rPWPart)
						w.resultChan <- outputZMWLine

					}
				}
			}
		}
	}
}
