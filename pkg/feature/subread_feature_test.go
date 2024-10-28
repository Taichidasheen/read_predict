package feature

import (
	"testing"
)

func Test_getCCS(t *testing.T) {

	seq1 := []byte{'A', 'A', 'C', 'C'}
	seq2 := []byte{'A', 'A', 'C', 'C'}
	seq3 := []byte{'A', 'A', 'C', 'A'}
	seq4 := []byte{'A', 'C', 'C', 'A'}
	seqs := [][]byte{seq1, seq2, seq3, seq4}
	cseq := getCCS(seqs)
	t.Logf("cseq:%s", string(cseq))
}

func Test_selectBaseByOrder(t *testing.T) {
	//maxCountBase := []byte{'T', 'C', 'G', 'A'}
	//maxCountBase := []byte{'C', 'T'}
	maxCountBase := []byte{'G', 'C', 'T', 'A'}

	base := selectBaseByOrder(maxCountBase)

	t.Logf("base:%v", base)

}
