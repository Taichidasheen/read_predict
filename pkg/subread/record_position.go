package subread

import (
	"fmt"
	"github.com/rs/zerolog/log"
	"strconv"
	"strings"
	"sync"
)

type RecordPositionDict struct {
	ChrMap               map[string]interface{}
	ChrCpgs              []string                      //chr1_cpg, chr2_cpg
	CpgLocationPositions map[string][]*LocatedPosition // chr1_cpg:[LocationPosition]
}

// TopRecordPositionDict 可以处理并写出的cpg
type TopRecordPositionDict struct {
	ChrCpgs              []string
	CpgLocationPositions map[string][]*LocatedPosition // chr1_cpg:[LocationPosition]
}

func processRecordPositions(recordPositionDict *RecordPositionDict, recordPositions chan *RecordPosition, refOrderMap map[string]int) (*TopRecordPositionDict, error) {
	//整理下上一轮的结果
	for recordPos := range recordPositions {
		chr := recordPos.RefChr
		if _, exist := refOrderMap[chr]; !exist {
			log.Warn().Msgf("unknown chr:%s, just skip", chr)
			continue
		}
		recordPositionDict.ChrMap[chr] = struct{}{}
		locatedCpgs := recordPos.LocatedCpgs
		for _, cpg := range locatedCpgs {
			chrCpg := fmt.Sprintf("%s_%d", chr, cpg)
			//log.Debug().Msgf("count:%d, chrCpg:%s", count, chrCpg)

			recordPositionDict.CpgLocationPositions[chrCpg] = append(recordPositionDict.CpgLocationPositions[chrCpg], recordPos.CpgLocatedPosition[cpg]) //dict.CpgLocationPositions
		}
	}
	// 对chrCpgs进行排序
	var chrCpgs []string
	for chrCpg, _ := range recordPositionDict.CpgLocationPositions {
		chrCpgs = append(chrCpgs, chrCpg)
	}
	sortedChrCpgs, err := sortChrCpgs(chrCpgs, refOrderMap)
	if err != nil {
		log.Error().Msgf("sortChrCpgs err:%v", err)
		return nil, err
	}
	recordPositionDict.ChrCpgs = sortedChrCpgs

	log.Debug().Msgf("len(sortedChrCpgs):%d, sortedChrCpgs:%v", len(sortedChrCpgs), sortedChrCpgs)

	//判断是不是可以安全的进行处理了
	if len(recordPositionDict.ChrMap) == 1 {
		//检查首尾距离
		headChrCpg := sortedChrCpgs[0]
		tailChrCpg := sortedChrCpgs[len(sortedChrCpgs)-1]

		//注意：这里没有检查split后的结果和err，因为前面sortChrCpgs方法已经检查过了
		headPos, _ := strconv.Atoi(strings.Split(headChrCpg, "_")[1])
		tailPos, _ := strconv.Atoi(strings.Split(tailChrCpg, "_")[1])

		log.Debug().Msgf("headChrCpg:%s, tailChrCpg:%s, headPos:%d, tailPos:%d", headChrCpg, tailChrCpg, headPos, tailPos)

		if tailPos-headPos > 50000 {
			//output the top 50% of the CpG in Feature_dict
			halfPos := len(sortedChrCpgs) / 2

			halfcpgs := sortedChrCpgs[0:halfPos]

			restCount := len(sortedChrCpgs) - len(halfcpgs)
			restCpgs := make([]string, restCount)

			numCopied := copy(restCpgs, sortedChrCpgs[halfPos:])
			if numCopied != restCount {
				log.Error().Msgf("array copy wrong, restCount:%d, numCopied:%d", restCount, numCopied)
				return nil, fmt.Errorf("array copy wrong")
			}

			log.Debug().Msgf("before delete len(recordPositionDict.ChrCpgs):%d", len(recordPositionDict.ChrCpgs))

			//记得从recordPositionDict中删除写出的部分
			recordPositionDict.ChrCpgs = restCpgs

			log.Debug().Msgf("after delete len(recordPositionDict.ChrCpgs):%d", len(recordPositionDict.ChrCpgs))

			cpgLocatedPositions := make(map[string][]*LocatedPosition)
			for _, chrCpg := range halfcpgs {
				cpgLocatedPositions[chrCpg] = recordPositionDict.CpgLocationPositions[chrCpg]
				//记得从recordPositionDict中删除写出的部分
				delete(recordPositionDict.CpgLocationPositions, chrCpg)
			}

			log.Debug().Msgf("after delete len(recordPositionDict.CpgLocationPositions):%d", len(recordPositionDict.CpgLocationPositions))

			//需要写出的部分
			topDict := &TopRecordPositionDict{
				ChrCpgs: halfcpgs,
			}
			topDict.CpgLocationPositions = cpgLocatedPositions

			return topDict, nil

		}
	}
	if len(recordPositionDict.ChrMap) > 1 {

		lastChrCpg := sortedChrCpgs[len(sortedChrCpgs)-1]
		lastChr := strings.Split(lastChrCpg, "_")[0]

		var frontChrCpgs []string
		cpgLocatedPositions := make(map[string][]*LocatedPosition)
		outputChr := make(map[string]interface{}) //记录下写出的chr，用于删除recordPositionDict.ChrMap
		for _, chrCpg := range sortedChrCpgs {
			chr := strings.Split(chrCpg, "_")[0]
			if chr == lastChr {
				break
			}
			outputChr[chr] = struct{}{}
			frontChrCpgs = append(frontChrCpgs, chrCpg)
			cpgLocatedPositions[chrCpg] = recordPositionDict.CpgLocationPositions[chrCpg]
			//记得从recordPositionDict中删除写出的部分
			delete(recordPositionDict.CpgLocationPositions, chrCpg)
		}

		//记得从recordPositionDict中删除写出的部分
		for key, _ := range outputChr {
			delete(recordPositionDict.ChrMap, key)
		}
		restCount := len(sortedChrCpgs) - len(frontChrCpgs)
		restCpgs := make([]string, restCount)

		numCopied := copy(restCpgs, sortedChrCpgs[len(frontChrCpgs):])
		if numCopied != restCount {
			log.Error().Msgf("array copy wrong, restCount:%d, numCopied:%d", restCount, numCopied)
			return nil, fmt.Errorf("array copy wrong")
		}
		//记得从recordPositionDict中删除写出的部分
		recordPositionDict.ChrCpgs = restCpgs

		//需要写出的部分
		topDict := &TopRecordPositionDict{
			ChrCpgs:              frontChrCpgs,
			CpgLocationPositions: cpgLocatedPositions,
		}
		return topDict, nil
	}

	log.Debug().Msgf("no data can output now, ChrCpgs:%v", len(recordPositionDict.ChrCpgs))

	return nil, nil

}

func processFinalRecordPositions(recordPositionDict *RecordPositionDict, recordPositions chan *RecordPosition, refOrderMap map[string]int) (*TopRecordPositionDict, error) {
	//整理下上一轮的结果
	for recordPos := range recordPositions {
		chr := recordPos.RefChr
		if _, exist := refOrderMap[chr]; !exist {
			log.Error().Msgf("unknown chr:%s, just skip", chr)
			continue
		}
		recordPositionDict.ChrMap[chr] = struct{}{}
		locatedCpgs := recordPos.LocatedCpgs

		for _, cpg := range locatedCpgs {
			chrCpg := fmt.Sprintf("%s_%d", chr, cpg)
			recordPositionDict.CpgLocationPositions[chrCpg] = append(recordPositionDict.CpgLocationPositions[chrCpg], recordPos.CpgLocatedPosition[cpg]) //dict.CpgLocationPositions
		}
	}
	// 对chrCpgs进行排序
	var chrCpgs []string
	for chrCpg, _ := range recordPositionDict.CpgLocationPositions {
		chrCpgs = append(chrCpgs, chrCpg)
	}
	sortedChrCpgs, err := sortChrCpgs(chrCpgs, refOrderMap)
	if err != nil {
		log.Error().Msgf("sortChrCpgs err:%v", err)
		return nil, err
	}
	recordPositionDict.ChrCpgs = sortedChrCpgs

	log.Debug().Msgf("processFinalRecordPositions sortedChrCpgs:%+v", sortedChrCpgs)

	//不用再判断是否可以安全处理，需要全部输出

	//不用再从recordPositionDict中删除写出的部分

	topDict := &TopRecordPositionDict{
		ChrCpgs:              recordPositionDict.ChrCpgs,
		CpgLocationPositions: recordPositionDict.CpgLocationPositions,
	}

	return topDict, nil

}

func sortCpgOutput2Chan(wg *sync.WaitGroup, cpgOutputChan chan *CpgOutput, resultChan chan string, refOrderMap map[string]int) {
	defer wg.Done()

	var cpgOutputs []*CpgOutput

	for cpgOutput := range cpgOutputChan {
		cpgOutputs = append(cpgOutputs, cpgOutput)
	}
	sortCpgOutput(cpgOutputs, refOrderMap)

	for _, cpgOutput := range cpgOutputs {
		for _, zmwLine := range cpgOutput.ZMWLines {
			resultChan <- zmwLine
		}
	}
}
