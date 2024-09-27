package subread

import (
	"fmt"
	"github.com/Taichidasheen/read_predict/pkg/feature"
	"github.com/Taichidasheen/read_predict/pkg/opt"
	"github.com/Taichidasheen/read_predict/pkg/record_flag"
	"github.com/Taichidasheen/read_predict/pkg/record_tag"
	"github.com/Taichidasheen/read_predict/pkg/util"
	"github.com/rs/zerolog/log"
	"strings"
)

type AlignedSubreadsFeatureWorker struct {
	locatedPositions []*LocatedPosition
	cpgOutputChan    chan *CpgOutput
	opts             opt.Options
	err              error
}

func NewAlignedSubreadsFeatureWorker(locatedPositions []*LocatedPosition, cpgOutputChan chan *CpgOutput, opts opt.Options) AlignedSubreadsFeatureWorker {
	return AlignedSubreadsFeatureWorker{
		locatedPositions: locatedPositions,
		cpgOutputChan:    cpgOutputChan,
		opts:             opts,
	}
}

func (w *AlignedSubreadsFeatureWorker) Task(num int) {
	cpgLocatedPositions := w.locatedPositions
	radius := w.opts.Radius
	scaleFlag := w.opts.ScaleFlag

	refChr := cpgLocatedPositions[0].RefChr
	locatedCpg := cpgLocatedPositions[0].LocatedCpg

	zmwLocatedPositions := make(map[string][]*LocatedPosition)
	//按照zmwname区分开
	for _, locatedPosition := range cpgLocatedPositions {
		readName := locatedPosition.Record.Name
		parts := strings.Split(readName, "/")
		zmwname := parts[0] + "_" + parts[1]
		zmwLocatedPositions[zmwname] = append(zmwLocatedPositions[zmwname], locatedPosition)
	}

	cpgOutput := &CpgOutput{
		Ref: refChr,
		Cpg: locatedCpg,
	}

	for zmwname, locatedPositions := range zmwLocatedPositions {

		zmwFeature := &feature.ZMWFeature{}

		zmwFeature.ZMWName = zmwname

		for _, locatedPos := range locatedPositions {

			record := locatedPos.Record
			recordTag, err := record_tag.ExtractSubreadsRecordTag(record)
			if err != nil {
				log.Warn().Msgf("extractSubreadsRecordTag err:%v", err)
				continue
			}
			readIpdList := recordTag.Ip
			readPwList := recordTag.Pw
			readSeqList := record.Seq.Expand()
			readIsReverse := record_flag.IsReverse(record.Flags)
			readQueryLength := len(readSeqList)
			readName := record.Name

			posOnSeq := locatedPos.PosOnSeq

			subreadFeature, err := feature.SubRead_cpg_K_Feature(posOnSeq, readIsReverse, radius, readQueryLength, readSeqList, readIpdList, readPwList)
			if err != nil {
				log.Error().Msgf("SubRead_cpg_K_Feature err:%v, readName:%s, locatedPos:%v", err, readName, locatedPos)
				continue
			}
			readQueryFeature := &feature.ReadQueryFeature{
				ReadQueryName:  readName,
				SubReadFeature: subreadFeature,
			}
			if subreadFeature.StrandSign == 'F' {
				zmwFeature.FStrandSubreadsFeature = append(zmwFeature.FStrandSubreadsFeature, readQueryFeature)
			}
			if subreadFeature.StrandSign == 'R' {
				zmwFeature.RStrandSubreadsFeature = append(zmwFeature.RStrandSubreadsFeature, readQueryFeature)
			}
		}

		zmwLine := formatZMWLine(zmwFeature, refChr, locatedCpg, scaleFlag)
		cpgOutput.ZMWLines = append(cpgOutput.ZMWLines, zmwLine)
	}

	log.Debug().Msgf("cpgOutput ref chr:%s, cpg:%d, len(zmwlines):%d", cpgOutput.Ref, cpgOutput.Cpg, len(cpgOutput.ZMWLines))

	//输出结果
	w.cpgOutputChan <- cpgOutput

}

type CpgOutput struct {
	Ref      string
	Cpg      int
	ZMWLines []string
}

// formatZMWLine 格式话每一行输出数据
func formatZMWLine(zmwFeature *feature.ZMWFeature, chr string, pos int, scaleFlag bool) string {

	zmwFSubDep := len(zmwFeature.FStrandSubreadsFeature)
	zmwRSubDep := len(zmwFeature.RStrandSubreadsFeature)

	zmwCCSFeature := feature.SumSubreadsZMW(zmwFeature, scaleFlag)

	fTemplateSeqStr, fIPDPart, fPWPart, rTemplateSeqStr, rIPDPart, rPWPart := "NA", "NA", "NA", "NA", "NA", "NA"

	if len(zmwCCSFeature.FCCSSeq) > 0 {
		fTemplateSeqStr = string(zmwCCSFeature.FCCSSeq)
	}
	if len(zmwCCSFeature.FCCSIPDList) > 0 {
		fIPDPart = util.FormatSlice(zmwCCSFeature.FCCSIPDList)
	}
	if len(zmwCCSFeature.FCCSPWList) > 0 {
		fPWPart = util.FormatSlice(zmwCCSFeature.FCCSPWList)
	}

	if len(zmwCCSFeature.RCCSSeq) > 0 {
		rTemplateSeqStr = string(zmwCCSFeature.RCCSSeq)
	}
	if len(zmwCCSFeature.RCCSIPDList) > 0 {
		rIPDPart = util.FormatSlice(zmwCCSFeature.RCCSIPDList)
	}
	if len(zmwCCSFeature.RCCSPWList) > 0 {
		rPWPart = util.FormatSlice(zmwCCSFeature.RCCSPWList)
	}

	zmwName := zmwFeature.ZMWName

	outputZMWLine := fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s",
		chr, pos, zmwName, zmwFSubDep, zmwRSubDep, fTemplateSeqStr, fIPDPart, fPWPart,
		rTemplateSeqStr, rIPDPart, rPWPart)
	return outputZMWLine
}
