package worker

import "github.com/biogo/hts/sam"

type AlignedSubreadsFeatureWorker struct {
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

func NewAlignedSubreadsFeatureWorker(record *sam.Record, cgListMap map[string][]int, resultChan chan string,
	mappingQ, minSubDep, maxSubDep, radius int, scaleFlag bool) AlignedSubreadsFeatureWorker {
	return AlignedSubreadsFeatureWorker{
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

func (w *AlignedSubreadsFeatureWorker) Task(num int) {

}
