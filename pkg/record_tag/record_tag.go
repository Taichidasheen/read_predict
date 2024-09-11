package record_tag

import (
	"fmt"
	"github.com/biogo/hts/sam"
	"log"
	"strings"
)

// RecordTag Notice:fn,rn这些tag在bam文件中是fn,rn(区分大小写)
type RecordTag struct {
	Fn int32
	Rn int32
	Fi []uint8
	Fp []uint8
	Ri []uint8
	Rp []uint8
	// HP:i:2   PS:i:11457
	HP string
	PS string

	Ip []uint8
	Pw []uint8
}

func ExtractRecordTag(record *sam.Record) (*RecordTag, error) {
	//fn
	fnAux, exist := record.Tag([]byte{'f', 'n'})
	if !exist {
		log.Printf("record fn not exist, record name:%s", record.Name)
		return nil, fmt.Errorf("fn not exist")
	}
	fn, ok := fnAux.Value().(int32)
	if !ok {
		log.Printf("record fn invalid, record name:%s", record.Name)
		return nil, fmt.Errorf("fn invalid")
	}

	//rn
	rnAux, exist := record.Tag([]byte{'r', 'n'})
	if !exist {
		log.Printf("record rn not exist, record name:%s", record.Name)
		return nil, fmt.Errorf("rn not exist")
	}
	rn, ok := rnAux.Value().(int32)
	if !ok {
		log.Printf("record rn invalid, record name:%s", record.Name)
		return nil, fmt.Errorf("rn invalid")
	}
	fmt.Println("record name:", record.Name, "fn:", fn, " rn:", rn)

	//fi
	fiAux, exist := record.Tag([]byte{'f', 'i'})
	if !exist {
		log.Printf("record fi not exist, record name:%s", record.Name)
		return nil, fmt.Errorf("fi not exist")
	}
	fi, ok := fiAux.Value().([]uint8)
	if !ok {
		log.Printf("record fi invalid, record name:%s", record.Name)
		return nil, fmt.Errorf("fi invalid")
	}

	//fp
	fpAux, exist := record.Tag([]byte{'f', 'p'})
	if !exist {
		log.Printf("record fp not exist, record name:%s", record.Name)
		return nil, fmt.Errorf("fp not exist")
	}
	fp, ok := fpAux.Value().([]uint8)
	if !ok {
		log.Printf("record fp invalid, record name:%s", record.Name)
		return nil, fmt.Errorf("fp invalid")
	}

	//ri
	riAux, exist := record.Tag([]byte{'r', 'i'})
	if !exist {
		log.Printf("record ri not exist, record name:%s", record.Name)
		return nil, fmt.Errorf("ri not exist")
	}
	ri, ok := riAux.Value().([]uint8)
	if !ok {
		log.Printf("record ri invalid, record name:%s", record.Name)
		return nil, fmt.Errorf("ri invalid")
	}

	//rp
	rpAux, exist := record.Tag([]byte{'r', 'p'})
	if !exist {
		log.Printf("record rp not exist, record name:%s", record.Name)
		return nil, fmt.Errorf("rp not exist")
	}
	rp, ok := rpAux.Value().([]uint8)
	if !ok {
		log.Printf("record rp invalid, record name:%s", record.Name)
		return nil, fmt.Errorf("rp invalid")
	}

	// HP:i:2   PS:i:11457
	var HP, PS string
	HPAux, exist := record.Tag([]byte{'H', 'P'})
	if !exist {
		log.Printf("record HP not exist, record name:%s", record.Name)
		HP = "X"
	} else {
		HPVal, ok := HPAux.Value().(int32)
		if !ok {
			log.Printf("record HP invalid, record name:%s", record.Name)
			return nil, fmt.Errorf("HP invalid")
		}
		HP = fmt.Sprintf("%d", HPVal)
	}

	PSAux, exist := record.Tag([]byte{'P', 'S'})
	if !exist {
		log.Printf("record PS not exist, record name:%s", record.Name)
		PS = "X"
	} else {
		PSVal, ok := PSAux.Value().(int32)
		if !ok {
			log.Printf("record PS invalid, record name:%s", record.Name)
			return nil, fmt.Errorf("PS invalid")
		}
		PS = fmt.Sprintf("%d", PSVal)
	}

	recordTag := &RecordTag{
		Fn: fn,
		Rn: rn,
		Fi: fi,
		Fp: fp,
		Ri: ri,
		Rp: rp,
		HP: HP,
		PS: PS,
	}
	return recordTag, nil
}

func ExtractSubreadsRecordTag(record *sam.Record) (*RecordTag, error) {
	//ip
	ipAux, exist := record.Tag([]byte{'i', 'p'})
	if !exist {
		log.Printf("record ip not exist, record name:%s", record.Name)
		return nil, fmt.Errorf("ip not exist")
	}
	ip, ok := ipAux.Value().([]uint8)
	if !ok {
		log.Printf("record ip invalid, record name:%s", record.Name)
		return nil, fmt.Errorf("ip invalid")
	}

	//pw
	pwAux, exist := record.Tag([]byte{'p', 'w'})
	if !exist {
		log.Printf("record pw not exist, record name:%s", record.Name)
		return nil, fmt.Errorf("pw not exist")
	}
	pw, ok := pwAux.Value().([]uint8)
	if !ok {
		log.Printf("record pw invalid, record name:%s", record.Name)
		return nil, fmt.Errorf("pw invalid")
	}

	recordTag := &RecordTag{
		Ip: ip,
		Pw: pw,
	}
	return recordTag, nil
}

func GenMLTag(probes []float32) (sam.Aux, error) {
	arrLen := len(probes)
	mlVal := make([]uint8, arrLen)
	for i, prob := range probes {
		pint := uint8(prob * 256)
		mlVal[i] = pint
	}
	return sam.NewAux(sam.NewTag("ML"), mlVal)
}

func GenMMTag(readSeqString string) (sam.Aux, error) {
	countC := 0
	var mmArr []string
	//注意：这里是len(readSeqString)-1,如果是len(readSeqString)可能会数组越界
	for i := 0; i < len(readSeqString)-1; i++ {
		if readSeqString[i] == 'C' && readSeqString[i:i+2] != "CG" {
			countC++
		}
		if readSeqString[i:i+2] == "CG" {
			mmArr = append(mmArr, fmt.Sprintf("%d", countC))
			countC = 0
		}
	}
	mmVal := "C+m," + strings.Join(mmArr, ",") + ";"
	return sam.NewAux(sam.NewTag("MM"), mmVal)
}

func GenAlignedMMTag(topStrand []byte, featurePosOnRead []int) (sam.Aux, error) {
	var mmArr []string

	//注意：这里是len(featurePos)-1,如果是len(featurePos)可能会数组越界
	for i := 0; i < len(featurePosOnRead)-1; i++ {
		start := featurePosOnRead[i]
		end := featurePosOnRead[i+1]
		countC := 0

		for j := start; j < end; j++ {
			if topStrand[j] == 'C' {
				countC++
			}
		}
		mmArr = append(mmArr, fmt.Sprintf("%d", countC-1))
	}
	mmVal := "C+m," + strings.Join(mmArr, ",") + ";"
	return sam.NewAux(sam.NewTag("MM"), mmVal)
}

func RemoveRecordTag(record *sam.Record) {
	auxes := record.AuxFields
	var newAuxes sam.AuxFields
	for _, aux := range auxes {
		tag := aux.Tag()
		if tag[0] == 'f' && tag[1] == 'i' {
			continue
		}
		if tag[0] == 'f' && tag[1] == 'p' {
			continue
		}
		if tag[0] == 'r' && tag[1] == 'i' {
			continue
		}
		if tag[0] == 'r' && tag[1] == 'p' {
			continue
		}
		newAuxes = append(newAuxes, aux)
	}
	record.AuxFields = newAuxes
}

func RemoveMMMLTag(record *sam.Record) {
	auxes := record.AuxFields
	var newAuxes sam.AuxFields
	for _, aux := range auxes {
		tag := aux.Tag()
		if tag[0] == 'M' && tag[1] == 'M' {
			continue
		}
		if tag[0] == 'M' && tag[1] == 'L' {
			continue
		}
		newAuxes = append(newAuxes, aux)
	}
	record.AuxFields = newAuxes
}
