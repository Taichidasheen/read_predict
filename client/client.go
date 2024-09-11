package main

import "fmt"

func main() {
	// 替换为你的FASTA文件路径
	filename := "/storage/yangjianLab/caoyujie/resource/human_ref/G_assemblies/GRCh38/v2_full/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa"

	records, err := ReadFASTA(filename)
	if err != nil {
		fmt.Println("Error reading file:", err)
		return
	}

	// 打印FASTA记录并查找 "CG" 出现的位置
	for _, record := range records {
		fmt.Println("Chromosome:", record.Chromosome)

		// 查找 "CG" 出现的位置
		cgPositions := record.FindCGPositions()
		fmt.Println("CG Positions:", cgPositions[:50])
	}
}
