package worker

import (
	"fmt"
	"github.com/Taichidasheen/read_predict/pkg/cpgpos"
	"github.com/Taichidasheen/read_predict/pkg/record_flag"
	"github.com/Taichidasheen/read_predict/pkg/record_tag"
	"github.com/biogo/hts/sam"
	"github.com/montanaflynn/stats"
	"github.com/rs/zerolog/log"
	"strings"
)

type LocateWorker struct {
	record     *sam.Record
	cgListMap  map[string][]int
	resultChan chan string
	mappingQ   int
	minSubDep  int
	maxSubDep  int
	radius     int
	scaleFlag  bool
	err        error
}

func NewLocateWorker(record *sam.Record, cgListMap map[string][]int, resultChan chan string,
	mappingQ, minSubDep, maxSubDep, radius int, scaleFlag bool) LocateWorker {
	return LocateWorker{
		record:     record,
		cgListMap:  cgListMap,
		resultChan: resultChan,
		mappingQ:   mappingQ,
		minSubDep:  minSubDep,
		maxSubDep:  maxSubDep,
		radius:     radius,
		scaleFlag:  scaleFlag,
	}
}

func (w *LocateWorker) Task(num int) {
	record := w.record
	radius := w.radius
	winsize := 2*radius + 1
	scaleFlag := w.scaleFlag
	if !record_flag.IsSecondary(record.Flags) && !record_flag.IsSupplementary(record.Flags) && int(record.MapQ) > w.mappingQ && record_flag.MatchingRatio(record) >= 0.1 {
		recordTag, err := record_tag.ExtractRecordTag(record)
		if err != nil {
			log.Error().Msgf("extractRecordTag err:%v", err)
			w.err = err
			return
		}
		fn := recordTag.Fn
		rn := recordTag.Rn
		totalSubreadsDep := fn + rn
		if totalSubreadsDep >= int32(w.minSubDep) && totalSubreadsDep <= int32(w.maxSubDep) && fn >= 1 && rn >= 1 {
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
				parts := strings.Split(readName, "/")
				zmwname := parts[0] + "_" + parts[1]
				haplotype := recordTag.HP
				haploTypeBlock := recordTag.PS

				locatedCpgs, cpgPosOnRead := cpgpos.LocateCpgPosOnSeq(alnRefStart, readCigar, overlappingCpg)
				log.Debug().Msgf("readName:%s, len(overlappingCpg):%d, len(locatedCpgs):%d", readName, len(overlappingCpg), len(locatedCpgs))
				log.Debug().Msgf("locatedCpgs:%+v", locatedCpgs)
				if len(locatedCpgs) >= 1 {
					for _, cpg := range locatedCpgs {
						posOnRead := cpgPosOnRead[cpg]
						//log.Debug().Msgf("readName:%s, cpg:%d, posOnRead:%d", readName, cpg, posOnRead)
						var fIPDList, fPWList, rIPDList, rPWList []float32
						var seqList []byte

						//this hifi read aligned to Ref_R_strand
						if readIsReverse {
							//(1)
							//For subreads mapped to Ref_R_strand, which carrying information when syn the complementary F_strand, key 'F' #
							//information is stored in fi, fp tags
							refFCOnFList := readQueryLength - posOnRead - 1
							leftRefFCOnFList := refFCOnFList - radius - 1
							rightRefFCOnFList := refFCOnFList + radius
							if leftRefFCOnFList > 0 && rightRefFCOnFList < readQueryLength {
								fIPDList, err = getkineticswin(readFiList[leftRefFCOnFList:rightRefFCOnFList], scaleFlag)
								if err != nil {
									log.Debug().Msgf("getkineticswin err:%+v", err)
									continue
								}
								fPWList, err = getkineticswin(readFpList[leftRefFCOnFList:rightRefFCOnFList], scaleFlag)
								if err != nil {
									log.Warn().Msgf("getkineticswin err:%+v", err)
									continue
								}
							}
							//(2)
							//For subreads mapped to Ref_F_strand, which carrying information when syn the complementary R_strand, key 'R' #
							//information is stored in ri, rp tags
							refRCOnRList := posOnRead + 1
							leftRefRCOnRList := refRCOnRList - radius - 1
							rightRefRCOnRList := refRCOnRList + radius
							if leftRefRCOnRList > 0 && rightRefRCOnRList < readQueryLength {
								rIPDList, err = getkineticswin(readRiList[leftRefRCOnRList:rightRefRCOnRList], scaleFlag)
								if err != nil {
									log.Warn().Msgf("getkineticswin err:%+v", err)
									continue
								}
								rPWList, err = getkineticswin(readRpList[leftRefRCOnRList:rightRefRCOnRList], scaleFlag)
								if err != nil {
									log.Warn().Msgf("getkineticswin err:%+v", err)
									continue
								}
								seqList = readSeqList[leftRefRCOnRList : rightRefRCOnRList+1]
								rIPDList = reverseSlice(rIPDList)
								rPWList = reverseSlice(rPWList)

							}

						} else {
							//(1)
							//For subreads mapped to Ref_F_strand, which carrying information when syn the complementary R_strand, key 'R' #
							//information is stored in fi, fp tags
							refRCOnFList := posOnRead + 1
							leftRefRCOnFList := refRCOnFList - radius - 1
							rightRefRCOnFList := refRCOnFList + radius
							if leftRefRCOnFList > 0 && rightRefRCOnFList < readQueryLength {
								rIPDList, err = getkineticswin(readFiList[leftRefRCOnFList:rightRefRCOnFList], scaleFlag)
								if err != nil {
									log.Warn().Msgf("getkineticswin err:%+v", err)
									continue
								}
								rPWList, err = getkineticswin(readFpList[leftRefRCOnFList:rightRefRCOnFList], scaleFlag)
								if err != nil {
									log.Warn().Msgf("getkineticswin err:%+v", err)
									continue
								}
								seqList = readSeqList[leftRefRCOnFList : rightRefRCOnFList+1]

							}
							//(2)
							//For subreads mapped to Ref_R_strand, which carrying information when syn the complementary F_strand, key 'F' #
							//information is stored in ri, rp tags
							refFCOnRList := readQueryLength - posOnRead - 1
							leftRefFCOnRList := refFCOnRList - radius - 1
							rightRefFCOnRList := refFCOnRList + radius
							if leftRefFCOnRList > 0 && rightRefFCOnRList < readQueryLength {
								fIPDList, err = getkineticswin(readRiList[leftRefFCOnRList:rightRefFCOnRList], scaleFlag)
								if err != nil {
									log.Warn().Msgf("getkineticswin err:%+v", err)
									continue
								}
								fPWList, err = getkineticswin(readRpList[leftRefFCOnRList:rightRefFCOnRList], scaleFlag)
								if err != nil {
									log.Warn().Msgf("getkineticswin err:%+v", err)
									continue
								}

								fIPDList = reverseSlice(fIPDList)
								fPWList = reverseSlice(fPWList)
							}
						}
						//输出
						ccsSeqString := string(seqList)
						fIPDPart := formatSlice(fIPDList)
						fPWPart := formatSlice(fPWList)
						rIPDPart := formatSlice(rIPDList)
						rPWPart := formatSlice(rPWList)

						if len(fIPDList) == winsize && len(fPWList) == winsize && len(rIPDList) == winsize && len(rPWList) == winsize {
							outputZMWLine := fmt.Sprintf("%s\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s",
								alnRefChr, cpg, ccsSeqString, zmwname, fn, rn, fIPDPart, fPWPart, rIPDPart, rPWPart, haplotype, haploTypeBlock)
							w.resultChan <- outputZMWLine

						}

					}
				}
			}
		}
	}
}

func getkineticswin(nums []uint8, scaleFlag bool) ([]float32, error) {

	if scaleFlag == false {
		var fnums []float32
		for _, num := range nums {
			fnums = append(fnums, float32(num))
		}
		return fnums, nil
	}
	//计算平均值和标准差
	if len(nums) == 0 {
		return nil, fmt.Errorf("empty input")
	}
	data := stats.LoadRawData(nums)
	mean, _ := stats.Mean(data)
	std, _ := stats.StandardDeviation(data)
	//fmt.Println("std:", std)
	//计算(x-mean)/std
	var kwin []float32
	for _, num := range nums {
		temp, err := stats.Round((float64(num)-mean)/std, 2)
		if err != nil {
			return nil, err
		}
		kwin = append(kwin, float32(temp))
	}
	return kwin, nil
}

func reverseSlice(arr []float32) []float32 {

	res := make([]float32, len(arr))
	length := len(arr)
	for i := length - 1; i >= 0; i-- {
		res[length-i-1] = arr[i]
	}
	return res
}

func reverseXReads(arr [][][]float32) [][][]float32 {

	res := make([][][]float32, len(arr))
	length := len(arr)
	for i := length - 1; i >= 0; i-- {
		res[length-i-1] = arr[i]
	}
	return res
}

func formatSlice[T any](arr []T) string {
	strs := make([]string, len(arr))
	for i, num := range arr {
		strs[i] = fmt.Sprintf("%v", num)
	}
	return strings.Join(strs, ",")
}
