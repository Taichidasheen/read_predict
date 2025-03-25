package feature

import (
	"fmt"
	"github.com/Taichidasheen/read_predict/pkg/util"
	"github.com/montanaflynn/stats"
	"log"
	"math"
)

type Feature struct {
	TemplateSeq     []byte
	TemplateIPDList []float32
	TemplatePWList  []float32

	ComTemplateSeq     []byte
	ComTemplateIPDList []float32
	ComTemplatePWList  []float32
}

func HiFiRead_cpg_K_Feature(posOnSeq int, readIsReverse bool, radius int, readQueryLength int, readSeqList []byte,
	readFiList, readFpList, readRiList, readRpList []uint8, scaleFlag bool) (*Feature, error) {
	/*start := time.Now()
	defer func() {
		log.Println("HiFiRead_cpg_K_Feature cost:", time.Since(start))
	}()*/
	/*
			2023-12-25:
			Related knowledge: DNA polymerase moves along in 3'->5' direction on the template, generating the 5'->3' reads.
			Return: (1) sequence context from the template of Ref_Forward_strand(should be in 3'->5', on reads is 5'->3') + IPD/PW on the reads while syning the template(should be in F_ipd_list and F_pw_list)
			            (2) sequence context from the template of Ref_Reverse_strand + IPD/PW on the reads while syning the template(should be in R_ipd_list and R_pw_list)
		    2024-01-09:
		    Related knowledge: DNA polymerase moves along in 3'->5' direction on the template, generating the 5'->3' reads.
		    Return: (1) TopHiFi_Template_seq in 3'->5' + IPD/PW while syning the TopHiFi_Template_seq
		                HiFiread_template_seq, HiFiread_template_ipd_list, HiFiread_template_pw_list
		            (2) ComHiFi_Template_seq in 3'->5' + IPD/PW while syning the ComHiFi_Template_seq
		                HiFiread_com_template_seq, HiFiread_com_template_ipd_list, HiFiread_com_template_pw_list
	*/

	//log.Printf("readName:%s, cpg:%d, posOnRead:%d", readName, cpg, posOnRead)
	var templateIPDList, templatePWList, comTemplateIPDList, comTemplatePWList []float32
	var templateSeq, comTemplateSeq []byte
	var err error

	//this hifi read aligned to Ref_R_strand
	if readIsReverse {
		//(1)
		//For subreads mapped to Ref_R_strand, which carrying information when syn the complementary F_strand, key 'F' #
		//information is stored in fi, fp tags
		refFCOnFList := readQueryLength - posOnSeq - 1
		//leftRefFCOnFList := refFCOnFList - radius - 1
		//rightRefFCOnFList := refFCOnFList + radius
		leftRefFCOnFList := refFCOnFList - radius
		rightRefFCOnFList := refFCOnFList + radius + 1
		if leftRefFCOnFList > 0 && rightRefFCOnFList < readQueryLength {
			templateIPDList, err = getkineticswinQuick(readFiList[leftRefFCOnFList:rightRefFCOnFList], scaleFlag)
			if err != nil {
				log.Printf("getkineticswin err:%+v", err)
				return nil, err
			}
			templatePWList, err = getkineticswinQuick(readFpList[leftRefFCOnFList:rightRefFCOnFList], scaleFlag)
			if err != nil {
				log.Printf("getkineticswin err:%+v", err)
				return nil, err
			}
			//fTemplateSeq = reverseSlice(readSeqList)[leftRefFCOnFList:rightRefFCOnFList]
			//fTemplateSeq = reverseSliceByte(readSeqList)[leftRefFCOnFList:rightRefFCOnFList]
			templateSeq = util.ReverseSliceByte(readSeqList[posOnSeq-radius : posOnSeq+radius+1])
		}
		//(2)
		//For subreads mapped to Ref_F_strand, which carrying information when syn the complementary R_strand, key 'R' #
		//information is stored in ri, rp tags
		refRCOnRList := posOnSeq + 1
		//leftRefRCOnRList := refRCOnRList - radius
		//rightRefRCOnRList := refRCOnRList + radius + 1
		leftRefRCOnRList := refRCOnRList - radius - 1
		rightRefRCOnRList := refRCOnRList + radius
		if leftRefRCOnRList > 0 && rightRefRCOnRList < readQueryLength {
			comTemplateIPDList, err = getkineticswinQuick(readRiList[leftRefRCOnRList:rightRefRCOnRList], scaleFlag)
			if err != nil {
				log.Printf("getkineticswin err:%+v", err)
				return nil, err
			}
			comTemplatePWList, err = getkineticswinQuick(readRpList[leftRefRCOnRList:rightRefRCOnRList], scaleFlag)
			if err != nil {
				log.Printf("getkineticswin err:%+v", err)
				return nil, err
			}
			comTemplateSeq = seqComplementary(readSeqList[leftRefRCOnRList+1 : rightRefRCOnRList+1])
			comTemplateIPDList = util.ReverseSlice(comTemplateIPDList)
			comTemplatePWList = util.ReverseSlice(comTemplatePWList)
		}

	} else { //this hifi read aligned to Ref_F_strand
		//(1)
		//For subreads mapped to Ref_F_strand, which carrying information when syn the complementary R_strand, key 'R' #
		//information is stored in fi, fp tags
		refRCOnFList := posOnSeq + 1
		leftRefRCOnFList := refRCOnFList - radius
		rightRefRCOnFList := refRCOnFList + radius + 1
		if leftRefRCOnFList > 0 && rightRefRCOnFList < readQueryLength {
			templateIPDList, err = getkineticswinQuick(readFiList[leftRefRCOnFList:rightRefRCOnFList], scaleFlag)
			if err != nil {
				log.Printf("getkineticswin err:%+v", err)
				return nil, err
			}
			templatePWList, err = getkineticswinQuick(readFpList[leftRefRCOnFList:rightRefRCOnFList], scaleFlag)
			if err != nil {
				log.Printf("getkineticswin err:%+v", err)
				return nil, err
			}
			templateSeq = seqComplementary(readSeqList[leftRefRCOnFList:rightRefRCOnFList])

		}
		//(2)
		//For subreads mapped to Ref_R_strand, which carrying information when syn the complementary F_strand, key 'F' #
		//information is stored in ri, rp tags
		refFCOnRList := readQueryLength - posOnSeq - 1
		leftRefFCOnRList := refFCOnRList - radius - 1
		rightRefFCOnRList := refFCOnRList + radius
		if leftRefFCOnRList > 0 && rightRefFCOnRList < readQueryLength {
			comTemplateIPDList, err = getkineticswinQuick(readRiList[leftRefFCOnRList:rightRefFCOnRList], scaleFlag)
			if err != nil {
				log.Printf("getkineticswin err:%+v", err)
				return nil, err
			}
			comTemplatePWList, err = getkineticswinQuick(readRpList[leftRefFCOnRList:rightRefFCOnRList], scaleFlag)
			if err != nil {
				log.Printf("getkineticswin err:%+v", err)
				return nil, err
			}
			comTemplateSeq = util.ReverseSliceByte(readSeqList[posOnSeq-radius : posOnSeq+radius+1])
			comTemplateIPDList = util.ReverseSlice(comTemplateIPDList)
			comTemplatePWList = util.ReverseSlice(comTemplatePWList)
		}
	}

	feature := &Feature{
		TemplateSeq:        templateSeq,
		TemplateIPDList:    templateIPDList,
		TemplatePWList:     templatePWList,
		ComTemplateSeq:     comTemplateSeq,
		ComTemplateIPDList: comTemplateIPDList,
		ComTemplatePWList:  comTemplatePWList,
	}
	if featureHasEmptyField(feature) {
		log.Printf("feature has empty field, posOnSeq:%d, feature:%+v", posOnSeq, feature)
		return nil, fmt.Errorf("empty feature field")
	}

	return feature, nil
}

func featureHasEmptyField(feature *Feature) bool {
	if len(feature.TemplateIPDList) == 0 || len(feature.TemplatePWList) == 0 || len(feature.TemplateSeq) == 0 ||
		len(feature.ComTemplateIPDList) == 0 || len(feature.ComTemplatePWList) == 0 || len(feature.ComTemplateSeq) == 0 {
		return true
	}
	return false
}

// 计算反向互补序列
func seqComplementary(seq []byte) []byte {
	complementedSeq := make([]byte, len(seq))
	for i := 0; i < len(seq); i++ {
		switch seq[i] {
		case 'A':
			complementedSeq[i] = 'T'
		case 'T':
			complementedSeq[i] = 'A'
		case 'C':
			complementedSeq[i] = 'G'
		case 'G':
			complementedSeq[i] = 'C'
		default:
			log.Println("invalid bytes:", seq[i])
		}
	}
	return complementedSeq
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

// getkineticswin的快速版实现
func getkineticswinQuick(nums []uint8, scaleFlag bool) ([]float32, error) {

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
		//log.Printf("std is zero, nums:%+v", nums)
		//return nil, fmt.Errorf("std is zero")
		//当std等于0时，设置成0.1
		std = 0.1
	}
	kwin := make([]float32, len(nums))
	//fmt.Println("std:", std)
	//计算(x-mean)/std
	for i, num := range nums {
		temp := round((float64(num)-mean)/std, 2)
		kwin[i] = float32(temp)
	}
	return kwin, nil
}

// Round 四舍五入，ROUND_HALF_UP 模式实现
// 返回将 val 根据指定精度 precision（十进制小数点后数字的数目）进行四舍五入的结果。precision 也可以是负数或零。
func round(val float64, precision int) float64 {
	if precision == 0 {
		return math.Round(val)
	}

	p := math.Pow10(precision)
	if precision < 0 {
		return math.Floor(val*p+0.5) * math.Pow10(-precision)
	}

	return math.Floor(val*p+0.5) / p
}

// 构造 13 x winsize矩阵
func formatClosedZMW(feature *Feature, npasses float32) [][]float32 {

	/*start := time.Now()
	defer func() {
		log.Println("formatClosedZMW cost:", time.Since(start))
	}()*/

	var matrix [][]float32
	fTemplateSeq := feature.TemplateSeq
	fIPDList := feature.TemplateIPDList
	fPWList := feature.TemplatePWList
	rTemplateSeq := feature.ComTemplateSeq
	rIPDList := feature.ComTemplateIPDList
	rPWList := feature.ComTemplatePWList

	winsize := len(fIPDList)
	if len(fTemplateSeq) != winsize || len(rTemplateSeq) != winsize {
		log.Printf("winsize:%d, fTemplateSeq:%s, rTemplateSeq:%s", winsize, string(fTemplateSeq), string(rTemplateSeq))
		return nil
	}
	seqDict := map[byte]int{'A': 0, 'T': 1, 'C': 2, 'G': 3}

	fseqMatArray := make([][]float32, 4)
	for i := 0; i < 4; i++ {
		fseqMatArray[i] = make([]float32, winsize)
	}

	rseqMatArray := make([][]float32, 4)
	for i := 0; i < 4; i++ {
		rseqMatArray[i] = make([]float32, winsize)
	}

	var npassesArr []float32

	//构造fseqMatArray，rseqMatArray，npassesArr
	for i := 0; i < winsize; i++ {
		fbase := fTemplateSeq[i]
		fseqMatArray[seqDict[fbase]][i] = 1
		rbase := rTemplateSeq[i]
		rseqMatArray[seqDict[rbase]][i] = 1
		npassesArr = append(npassesArr, npasses)
	}

	//13 x winsize matrix
	matrix = append(matrix, fseqMatArray...)
	matrix = append(matrix, fIPDList, fPWList)
	matrix = append(matrix, rseqMatArray...)
	matrix = append(matrix, rIPDList, rPWList, npassesArr)

	return matrix
}

func transpose2D(plane [][]float32) [][]float32 {
	rows := len(plane)
	cols := len(plane[0])

	result := make([][]float32, cols)
	for i := range result {
		result[i] = make([]float32, rows)
	}

	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			result[j][i] = plane[i][j]
		}
	}
	return result
}

// 实现 np.transpose(X_reads, (0, 2, 1))
func transpose3D(cube [][][]float32) [][][]float32 {
	/*start := time.Now()
	defer func() {
		log.Println("transpose3D cost:", time.Since(start))
	}()*/

	var result [][][]float32
	for _, plane := range cube {
		posedPlane := transpose2D(plane)
		result = append(result, posedPlane)
	}
	return result
}

// 构造 winsize x 13 矩阵
func FormatTransposedClosedZMW(feature *Feature, npasses float32) [][]float32 {

	/*start := time.Now()
	defer func() {
		log.Println("formatTransposedClosedZMW cost:", time.Since(start))
	}()*/

	fTemplateSeq := feature.TemplateSeq
	fIPDList := feature.TemplateIPDList
	fPWList := feature.TemplatePWList
	rTemplateSeq := feature.ComTemplateSeq
	rIPDList := feature.ComTemplateIPDList
	rPWList := feature.ComTemplatePWList

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
