package feature

import (
	"fmt"
	"github.com/Taichidasheen/read_predict/pkg/util"
	"log"
	"math"
)

type SubReadFeature struct {
	IPDList     []uint8
	PWList      []uint8
	TemplateSeq []byte
	StrandSign  byte
}

type ZMWFeature struct {
	ZMWName                string
	FStrandSubreadsFeature []*ReadQueryFeature
	RStrandSubreadsFeature []*ReadQueryFeature
}

type ReadQueryFeature struct {
	ReadQueryName string
	*SubReadFeature
}

type ZMWCCSFeature struct {
	FCCSSeq     []byte
	FCCSIPDList []float32
	FCCSPWList  []float32

	RCCSSeq     []byte
	RCCSIPDList []float32
	RCCSPWList  []float32
}

func subFeatureHasEmptyField(feature *SubReadFeature) bool {
	if len(feature.IPDList) == 0 || len(feature.PWList) == 0 || len(feature.TemplateSeq) == 0 {
		return true
	}
	return false
}

/*
	 	2024-01-03: subreads keep raw ipd/pw for each of them
		2024-01-02: return Seq_context_on_Template(3'->5'), ipd, pw
*/
func SubRead_cpg_K_Feature(posOnSeq int, readIsReverse bool, radius int, readQueryLength int, readSeqList []byte,
	readIPDList, readPWList []uint8) (*SubReadFeature, error) {
	var ipdList, pwList []uint8
	var templateSeq []byte
	var strandSign byte

	//# subreads mapped to Ref(-)_strand, synthesising the Ref(+)_template, ipd/pw are therefor with respected to Ref(+). Template is Ref(+), but should be reverse into 3'-5'
	if readIsReverse {
		//5'->3' position on read
		adjustPosOnRead := readQueryLength - posOnSeq - 1
		lPos := adjustPosOnRead - radius - 1
		rPos := adjustPosOnRead + radius
		strandSign = 'F' //'Ref(+)'#Template strand sign
		if lPos > 0 && rPos < readQueryLength {
			ipdList = readIPDList[lPos:rPos]
			pwList = readPWList[lPos:rPos]
			templateSeq = util.ReverseSliceByte(readSeqList[posOnSeq-radius : posOnSeq+radius+1])
		} else {
			log.Printf("F posOnSeq:%d, readQueryLength:%d, lPos:%d, rPos:%d", posOnSeq, readQueryLength, lPos, rPos)
		}

	} else {
		//subreads mapped to Ref(+), synthesising the Ref(-)_template, ipd/pw are therefor with respected to Ref(-). Template is Ref(-), and should be adjusted to 3'->5'
		//5'->3' position on read
		adjustPosOnRead := posOnSeq + 1
		lPos := adjustPosOnRead - radius
		rPos := adjustPosOnRead + radius + 1
		strandSign = 'R' //'Ref(-)'
		if lPos > 0 && rPos < readQueryLength {
			ipdList = readIPDList[lPos:rPos]
			pwList = readPWList[lPos:rPos]
			templateSeq = seqComplementary(readSeqList[lPos:rPos])
		} else {
			log.Printf("R posOnSeq:%d, readQueryLength:%d, lPos:%d, rPos:%d", posOnSeq, readQueryLength, lPos, rPos)
		}
	}
	sfeature := &SubReadFeature{
		TemplateSeq: templateSeq,
		IPDList:     ipdList,
		PWList:      pwList,
		StrandSign:  strandSign,
	}
	if subFeatureHasEmptyField(sfeature) {
		log.Printf("sub feature has empty field, posOnSeq:%d, sfeature:%+v", posOnSeq, sfeature)
		return nil, fmt.Errorf("empty feature field")
	}
	return sfeature, nil
}

// SumSubreadsZMW 生成ccs feature
func SumSubreadsZMW(zmwFeature *ZMWFeature, scaleFlag bool) *ZMWCCSFeature {

	var fseqs, rseqs [][]byte
	var fipds, ripds [][]uint8
	var fpws, rpws [][]uint8
	for _, feature := range zmwFeature.FStrandSubreadsFeature {
		fseqs = append(fseqs, feature.TemplateSeq)
		fipds = append(fipds, feature.IPDList)
		fpws = append(fpws, feature.PWList)
	}

	for _, feature := range zmwFeature.RStrandSubreadsFeature {
		rseqs = append(rseqs, feature.TemplateSeq)
		ripds = append(ripds, feature.IPDList)
		rpws = append(rpws, feature.PWList)
	}

	//ccs seq
	fccsSeq := getCCS(fseqs)
	rccsSeq := getCCS(rseqs)

	//average ipd, pw
	favgIPD := getColumnsAverage(fipds)
	favgPW := getColumnsAverage(fpws)

	ravgIPD := getColumnsAverage(ripds)
	ravgPW := getColumnsAverage(rpws)

	if scaleFlag {
		favgIPD = getScaleList(favgIPD)
		favgPW = getScaleList(favgPW)
		ravgIPD = getScaleList(ravgIPD)
		ravgPW = getScaleList(ravgPW)
	}

	zmwCCSFeature := &ZMWCCSFeature{
		FCCSSeq:     fccsSeq,
		FCCSIPDList: favgIPD,
		FCCSPWList:  favgPW,

		RCCSSeq:     rccsSeq,
		RCCSIPDList: ravgIPD,
		RCCSPWList:  ravgPW,
	}
	return zmwCCSFeature

}

func getCCS(seqs [][]byte) []byte {

	rows := len(seqs)
	if rows == 0 {
		return nil
	}

	numColumns := len(seqs[0])
	cSeq := make([]byte, numColumns)

	for col := 0; col < numColumns; col++ {
		counts := make(map[byte]int)

		for _, row := range seqs {
			counts[row[col]]++
		}

		var maxCount int
		var maxCountBase []byte

		//统计出出现次数最多的base
		for element, count := range counts {
			if count > maxCount {
				maxCount = count
				maxCountBase = []byte{element}
			}
			if count == maxCount {
				maxCount = count
				maxCountBase = append(maxCountBase, element)
			}
		}
		/*var cBase byte
		if len(maxCountBase) > 1 {
			//随机选择一个base
			idx := rand.Intn(len(maxCountBase))
			cBase = maxCountBase[idx]
		}
		if len(maxCountBase) == 1 {
			cBase = maxCountBase[0]
		}*/
		//去掉随机逻辑，按照ATCG的顺序进行选择
		cBase := selectBaseByOrder(maxCountBase)
		cSeq[col] = cBase
	}
	return cSeq
}

var baseOrderMap = map[byte]int{
	'A': 4,
	'T': 3,
	'C': 2,
	'G': 1,
}

func selectBaseByOrder(maxCountBase []byte) byte {
	selectedBase := maxCountBase[0]
	for _, base := range maxCountBase {
		if baseOrderMap[base] > baseOrderMap[selectedBase] {
			selectedBase = base
		}
	}
	return selectedBase
}

func getColumnsAverage(matrix [][]uint8) []float32 {
	rows := len(matrix)
	if rows == 0 {
		return nil
	}
	columns := len(matrix[0])
	colAvg := make([]float32, columns)

	for j := 0; j < columns; j++ {
		var sum int
		for i := 0; i < rows; i++ {
			sum = sum + int(matrix[i][j])
		}
		colAvg[j] = float32(sum / rows)
	}
	return colAvg
}

// getScaleList
func getScaleList(nums []float32) []float32 {
	//计算平均值和标准差
	if len(nums) == 0 {
		return nil
	}
	arrLen := float64(len(nums))
	var sum, mean, std float64
	for _, num := range nums {
		sum += float64(num)
	}
	mean = sum / arrLen
	for _, num := range nums {
		std += math.Pow(float64(num)-mean, 2)
	}
	std = math.Sqrt(std / arrLen)
	if std == 0 {
		return nums
	}
	kwin := make([]float32, len(nums))
	//fmt.Println("std:", std)
	//计算(x-mean)/std
	for i, num := range nums {
		temp := round((float64(num)-mean)/std, 2)
		kwin[i] = float32(temp)
	}
	return kwin
}

type CpgPredictFeature struct {
	ClosedZMWX        [][][]float32
	ClosedNames       []string
	ClosedNpassesList []int

	OpenZMWX        [][][]float32
	OpenNames       []string
	OpenNpassesList []int
}

type ZMWPredictFeature struct {
	Type    string
	Matrix  [][]float32
	Name    string
	Npasses int
}

func GenZMWPredictFeature(zmwFeature *ZMWFeature, scaleFlag bool) *ZMWPredictFeature {
	zmwCCSFeature := SumSubreadsZMW(zmwFeature, scaleFlag)
	zmwPredictFeature := &ZMWPredictFeature{}
	zmwPredictFeature.Name = zmwFeature.ZMWName

	lenFStrand := len(zmwFeature.FStrandSubreadsFeature)
	lenRStrand := len(zmwFeature.RStrandSubreadsFeature)

	if lenFStrand > 0 && lenRStrand > 0 {
		zmwPredictFeature.Type = "Closed"
		npasses := lenFStrand + lenRStrand
		matrix := formatTransposedClosedZMWCCS(zmwCCSFeature, float32(npasses))
		zmwPredictFeature.Matrix = matrix
		zmwPredictFeature.Npasses = npasses
	}
	if lenFStrand > 0 && lenRStrand == 0 {
		zmwPredictFeature.Type = "Open"
		npass := lenFStrand
		matrix := formatTransposedOpenZMWCCS(zmwCCSFeature.FCCSSeq, zmwCCSFeature.FCCSIPDList, zmwCCSFeature.FCCSPWList, float32(npass))
		zmwPredictFeature.Matrix = matrix
		zmwPredictFeature.Npasses = npass
	}
	if lenFStrand == 0 && lenRStrand > 0 {
		zmwPredictFeature.Type = "Open"
		npasses := lenRStrand
		matrix := formatTransposedOpenZMWCCS(zmwCCSFeature.RCCSSeq, zmwCCSFeature.RCCSIPDList, zmwCCSFeature.RCCSPWList, float32(npasses))
		zmwPredictFeature.Matrix = matrix
		zmwPredictFeature.Npasses = npasses
	}

	return zmwPredictFeature
}

// 构造 winsize x 13 矩阵
func formatTransposedClosedZMWCCS(zmwCCSFeature *ZMWCCSFeature, npasses float32) [][]float32 {

	/*start := time.Now()
	defer func() {
		log.Println("formatTransposedClosedZMW cost:", time.Since(start))
	}()*/

	fTemplateSeq := zmwCCSFeature.FCCSSeq
	fIPDList := zmwCCSFeature.FCCSIPDList
	fPWList := zmwCCSFeature.FCCSPWList
	rTemplateSeq := zmwCCSFeature.RCCSSeq
	rIPDList := zmwCCSFeature.RCCSIPDList
	rPWList := zmwCCSFeature.RCCSPWList

	winsize := len(fIPDList)
	if len(fTemplateSeq) != winsize || len(rTemplateSeq) != winsize {
		log.Printf("winsize:%d, fTemplateSeq:%s, rTemplateSeq:%s", winsize, string(fTemplateSeq), string(rTemplateSeq))
		return nil
	}
	seqDict := map[byte]int{'A': 0, 'T': 1, 'C': 2, 'G': 3}

	//winsize x 13 matrix
	matrix := make([][]float32, winsize)
	for i, _ := range matrix {
		matrix[i] = make([]float32, 13)
	}

	for i := 0; i < winsize; i++ {
		//ATCG, fi, fp
		fbase := fTemplateSeq[i]
		matrix[i][seqDict[fbase]] = 1
		matrix[i][4] = fIPDList[i]
		matrix[i][5] = fPWList[i]

		//ATCG, ri, rp
		rbase := rTemplateSeq[i]
		matrix[i][6+seqDict[rbase]] = 1
		matrix[i][10] = rIPDList[i]
		matrix[i][11] = rPWList[i]
		matrix[i][12] = npasses
	}

	return matrix
}

// 构造 winsize x 7 矩阵
func formatTransposedOpenZMWCCS(templateSeq []byte, ipdList []float32, pwList []float32, npasses float32) [][]float32 {

	/*start := time.Now()
	defer func() {
		log.Println("formatTransposedClosedZMW cost:", time.Since(start))
	}()*/

	winsize := len(ipdList)
	if len(templateSeq) != winsize {
		log.Printf("winsize:%d, templateSeq:%s", winsize, string(templateSeq))
		return nil
	}
	seqDict := map[byte]int{'A': 0, 'T': 1, 'C': 2, 'G': 3}

	//winsize x 7 matrix
	matrix := make([][]float32, winsize)
	for i, _ := range matrix {
		matrix[i] = make([]float32, 7)
	}

	for i := 0; i < winsize; i++ {
		//ATCG, fi, fp
		fbase := templateSeq[i]
		matrix[i][seqDict[fbase]] = 1
		matrix[i][4] = ipdList[i]
		matrix[i][5] = pwList[i]
		matrix[i][6] = npasses
	}

	return matrix
}
