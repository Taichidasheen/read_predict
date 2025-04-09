package subread

import (
	"github.com/Taichidasheen/read_predict/pkg/cpgpos"
	"github.com/Taichidasheen/read_predict/pkg/opt"
	"github.com/Taichidasheen/read_predict/pkg/record_flag"
	"github.com/biogo/hts/sam"
	"github.com/rs/zerolog/log"
)

// RecordPosition 存放一个record所有命中的cpg
type RecordPosition struct {
	RefChr             string
	LocatedCpgs        []int
	CpgLocatedPosition map[int]*LocatedPosition
}

// LocatedPosition 存放record单个命中的cpg的信息
type LocatedPosition struct {
	Record     *sam.Record
	RefChr     string
	LocatedCpg int
	PosOnSeq   int
}

func FindLocatedCpgs(record *sam.Record, opts opt.Options, cgListMap map[string][]int) *RecordPosition {

	mappingQ := opts.MappingQ
	radius := opts.Radius

	alnRefChr := record.Ref.Name()
	alnRefStart := record.Pos
	alnRefEnd := record.End()
	log.Debug().Msgf("readName:%s, alnRefStart:%d, alnRefEnd:%d", record.Name, alnRefStart, alnRefEnd)
	cgList := cgListMap[alnRefChr]

	if !record_flag.IsSecondary(record.Flags) && !record_flag.IsSupplementary(record.Flags) && int(record.MapQ) > mappingQ && record_flag.MatchingRatio(record) >= 0.1 {

		overlappingCpg := cpgpos.FindOverlappingCpg(cgList, alnRefStart, alnRefEnd)
		//log.Debug().Msgf("count:%d, overlappingCpg:%+v", count, overlappingCpg)

		if len(overlappingCpg) >= 1 {
			readCigar := record.Cigar
			readName := record.Name
			readSeqList := record.Seq.Expand()
			readQueryLength := len(readSeqList)

			locatedCpgs, cpgPosOnSeq := cpgpos.LocateCpgPosOnSeq(alnRefStart, readCigar, overlappingCpg)
			log.Debug().Msgf("readName:%s, len(overlappingCpg):%d, len(locatedCpgs):%d", readName, len(overlappingCpg), len(locatedCpgs))
			//log.Debug().Msgf("count:%d, locatedCpgs:%+v", count, locatedCpgs)

			if len(locatedCpgs) >= 1 {
				var filteredCpgs []int
				cpgLocatedPosition := make(map[int]*LocatedPosition)
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
					filteredCpgs = append(filteredCpgs, cpg)
					locPos := &LocatedPosition{
						Record:     record,
						RefChr:     alnRefChr,
						LocatedCpg: cpg,
						PosOnSeq:   cpgPosOnSeq[cpg],
					}
					cpgLocatedPosition[cpg] = locPos
				}
				if len(filteredCpgs) > 0 {
					recordPos := &RecordPosition{
						RefChr:             alnRefChr,
						LocatedCpgs:        filteredCpgs,
						CpgLocatedPosition: cpgLocatedPosition,
					}
					return recordPos
				}
			}
		}
	}
	return nil
}
