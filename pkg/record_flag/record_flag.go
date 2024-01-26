package record_flag

import "github.com/biogo/hts/sam"

func IsReverse(flag sam.Flags) bool {
	return flag&sam.Reverse == sam.Reverse
}

func IsSecondary(flag sam.Flags) bool {
	return flag&sam.Secondary == sam.Secondary
}

func IsSupplementary(flag sam.Flags) bool {
	return flag&sam.Supplementary == sam.Supplementary
}

func MatchingRatio(record *sam.Record) float32 {
	cigars := record.Cigar
	totalLen := 0
	matchingLen := 0
	for _, cigar := range cigars {
		op := cigar.Type()
		count := cigar.Len()
		// 0:M, 1:I, 4:S, 7:=, 8:X ||||||||||||| sum of "M/I/S/=/X" is the SEQ length
		if op == sam.CigarMatch || op == sam.CigarInsertion || op == sam.CigarSoftClipped ||
			op == sam.CigarEqual || op == sam.CigarMismatch {
			totalLen += count
		}
		if op == sam.CigarMatch || op == sam.CigarEqual || op == sam.CigarMismatch {
			matchingLen += count
		}
	}
	ratio := float32(matchingLen) / float32(totalLen)
	return ratio
}
