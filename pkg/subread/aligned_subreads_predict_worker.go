package subread

import (
	"fmt"
	"github.com/Taichidasheen/read_predict/pkg/feature"
	"github.com/Taichidasheen/read_predict/pkg/opt"
	"github.com/Taichidasheen/read_predict/pkg/predict"
	"github.com/Taichidasheen/read_predict/pkg/record_flag"
	"github.com/Taichidasheen/read_predict/pkg/record_tag"
	tf "github.com/wamuir/graft/tensorflow"
	"log"
	"strings"
)

type AlignedSubreadsPredictWorker struct {
	closedModel      *tf.SavedModel
	openModel        *tf.SavedModel
	locatedPositions []*LocatedPosition
	cpgOutputChan    chan *CpgOutput
	opts             opt.Options
	err              error
}

func NewAlignedSubreadsPredictWorker(closedModel, openModel *tf.SavedModel,
	locatedPositions []*LocatedPosition, cpgOutputChan chan *CpgOutput, opts opt.Options) AlignedSubreadsPredictWorker {
	return AlignedSubreadsPredictWorker{
		closedModel:      closedModel,
		openModel:        openModel,
		locatedPositions: locatedPositions,
		cpgOutputChan:    cpgOutputChan,
		opts:             opts,
	}
}

func (w *AlignedSubreadsPredictWorker) Task(num int) {
	closedModel := w.closedModel
	openModel := w.openModel
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

	cpgPredictFeature := &feature.CpgPredictFeature{}
	for zmwname, locatedPositions := range zmwLocatedPositions {

		zmwFeature := &feature.ZMWFeature{}

		zmwFeature.ZMWName = zmwname

		for _, locatedPos := range locatedPositions {

			record := locatedPos.Record
			recordTag, err := record_tag.ExtractSubreadsRecordTag(record)
			if err != nil {
				log.Printf("extractSubreadsRecordTag err:%v", err)
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
				log.Printf("SubRead_cpg_K_Feature err:%v, readName:%s, locatedPos:%v", err, readName, locatedPos)
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

		zmwPredictFeature := feature.GenZMWPredictFeature(zmwFeature, scaleFlag)
		if zmwPredictFeature.Type == "Closed" {
			cpgPredictFeature.ClosedZMWX = append(cpgPredictFeature.ClosedZMWX, zmwPredictFeature.Matrix)
			cpgPredictFeature.ClosedNames = append(cpgPredictFeature.ClosedNames, zmwPredictFeature.Name)
			cpgPredictFeature.ClosedNpassesList = append(cpgPredictFeature.ClosedNpassesList, zmwPredictFeature.Npasses)
		}
		if zmwPredictFeature.Type == "Open" {
			cpgPredictFeature.OpenZMWX = append(cpgPredictFeature.OpenZMWX, zmwPredictFeature.Matrix)
			cpgPredictFeature.OpenNames = append(cpgPredictFeature.OpenNames, zmwPredictFeature.Name)
			cpgPredictFeature.OpenNpassesList = append(cpgPredictFeature.OpenNpassesList, zmwPredictFeature.Npasses)
		}
	}

	//predict
	//closed predict
	if len(cpgPredictFeature.ClosedZMWX) > 0 {
		closedProbes, err := predict.Predict(closedModel, cpgPredictFeature.ClosedZMWX)
		if err != nil {
			log.Printf("closed predict err:%+v, input:%+v", err, cpgPredictFeature.ClosedZMWX)
			w.err = err
			return
		}
		//构造输出
		predictLine, err := formatZMWPredictLine(refChr, locatedCpg, cpgPredictFeature.ClosedNames, cpgPredictFeature.ClosedNpassesList, closedProbes, "C")
		if err != nil {
			log.Printf("formatZMWPredictLine err:%+v", err)
			w.err = err
			return
		}
		cpgOutput.ZMWLines = append(cpgOutput.ZMWLines, predictLine)
	}

	//open predict
	if len(cpgPredictFeature.OpenZMWX) > 0 {
		openProbes, err := predict.Predict(openModel, cpgPredictFeature.OpenZMWX)
		if err != nil {
			log.Printf("predict err:%+v, input:%+v", err, cpgPredictFeature.OpenZMWX)
			w.err = err
			return
		}
		//构造输出
		predictLine, err := formatZMWPredictLine(refChr, locatedCpg, cpgPredictFeature.OpenNames, cpgPredictFeature.OpenNpassesList, openProbes, "0")
		if err != nil {
			log.Printf("formatZMWPredictLine err:%+v", err)
			w.err = err
			return
		}
		cpgOutput.ZMWLines = append(cpgOutput.ZMWLines, predictLine)
	}

	log.Printf("cpgOutput ref chr:%s, cpg:%d, len(zmwlines):%d", cpgOutput.Ref, cpgOutput.Cpg, len(cpgOutput.ZMWLines))

	//输出结果
	w.cpgOutputChan <- cpgOutput

}

// formatZMWData ZMW_Type = 'C' or 'O'
func formatZMWPredictLine(chr string, pos int, zmwNames []string, npasses []int, probes []float32, zmwType string) (string, error) {

	if len(zmwNames) != len(probes) || len(npasses) != len(probes) {
		log.Printf("non uniform length, zmwNames:%v, npasses:%v, probes:%v", zmwNames, npasses, probes)
		return "", fmt.Errorf("non uniform length")
	}

	allLines := make([]string, len(probes))

	for i, prob := range probes {
		var bi int
		if prob >= 0.5 {
			bi = 1
		} else {
			bi = 0
		}
		zmwName := zmwNames[i]
		subDepth := npasses[i]
		hap := "X"
		hapBlk := "X"
		strand := "+/-"
		line := fmt.Sprintf("%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%v\t%d",
			chr, pos, zmwName, zmwType, subDepth, hap, hapBlk, strand, probes[i], bi)
		allLines[i] = line
	}
	allLinesStr := strings.Join(allLines, "\n")

	return allLinesStr, nil

}
