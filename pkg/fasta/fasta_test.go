package fasta

import (
	"fmt"
	"testing"
)

func TestReadFASTA(t *testing.T) {

	// 替换为你的FASTA文件路径
	filename := "/storage/yangjianLab/caoyujie/resource/human_ref/G_assemblies/GRCh38/v2_full/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa"

	records, err := ReadFASTA(filename)
	if err != nil {
		t.Log("Error reading file:", err)
		return
	}

	// 打印FASTA记录
	for _, record := range records {
		fmt.Println("Chromosome:", record.Chromosome)

		t.Log("Sequence:", record.Sequence[:50], "...") // 只显示前50个字符
	}

}
