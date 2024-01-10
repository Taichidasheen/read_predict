package worker

import (
	"testing"
)

func Test_transpose2D(t *testing.T) {
	reads := [][]float32{
		{1, 2, 3},
		{4, 5, 6},
		{7, 8, 9},
	}
	res := transpose2D(reads)
	t.Logf("res:%v", res)
}

func Test_formatClosedZMW(t *testing.T) {

	ftemplateseq := []byte("TTTGTT")
	fIPDList := []float32{-0.4, -0.48, -0.4, -0.33, 0.27, -0.18}
	fPWList := []float32{-0.48, -0.21, 0.06, -0.3, -0.03, -0.03}
	rtemplateseq := []byte("AACAAA")
	rIPDList := []float32{-0.3, -0.68, -0.65, -0.61, -0.72, -0.51}
	rPWList := []float32{-0.16, -0.23, -0.6, 0.21, -0.75, -0.75}
	npasses := float32(5)

	feature := &Feature{
		TemplateSeq:        ftemplateseq,
		TemplateIPDList:    fIPDList,
		TemplatePWList:     fPWList,
		ComTemplateSeq:     rtemplateseq,
		ComTemplateIPDList: rIPDList,
		ComTemplatePWList:  rPWList,
	}
	matrix := formatClosedZMW(feature, npasses)
	t.Logf("matrix:%+v", matrix)

	transpose := transpose2D(matrix)
	t.Logf("transpose matrix:%+v", transpose)

}

func Test_formatTransposedClosedZMW(t *testing.T) {
	ftemplateseq := []byte("TTTGTT")
	fIPDList := []float32{-0.4, -0.48, -0.4, -0.33, 0.27, -0.18}
	fPWList := []float32{-0.48, -0.21, 0.06, -0.3, -0.03, -0.03}
	rtemplateseq := []byte("AACAAA")
	rIPDList := []float32{-0.3, -0.68, -0.65, -0.61, -0.72, -0.51}
	rPWList := []float32{-0.16, -0.23, -0.6, 0.21, -0.75, -0.75}
	npasses := float32(5)

	feature := &Feature{
		TemplateSeq:        ftemplateseq,
		TemplateIPDList:    fIPDList,
		TemplatePWList:     fPWList,
		ComTemplateSeq:     rtemplateseq,
		ComTemplateIPDList: rIPDList,
		ComTemplatePWList:  rPWList,
	}
	transpose := formatTransposedClosedZMW(feature, npasses)
	t.Logf("transpose matrix:%+v", transpose)

}
