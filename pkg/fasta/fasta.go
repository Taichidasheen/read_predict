package fasta

import (
	"bufio"
	"os"
	"strings"
)

// FASTARecord 用于存储每个FASTA记录
type FASTARecord struct {
	Chromosome string // 染色体标识符（例如: "chr1"）
	Sequence   string // 生物序列
}

// ReadFASTA 读取FASTA文件并解析记录
func ReadFASTA(filename string) ([]*FASTARecord, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var records []*FASTARecord
	scanner := bufio.NewScanner(file)
	var currentRecord *FASTARecord
	var sequenceBuilder strings.Builder // 使用 strings.Builder 优化字符串拼接

	for scanner.Scan() {
		line := scanner.Text()

		//fmt.Println(line)

		// 判断是否是序列标识符行（以'>'开头）
		if strings.HasPrefix(line, ">") {
			if currentRecord != nil {
				// 将当前记录添加到结果集中
				currentRecord.Sequence = sequenceBuilder.String() // 获取拼接后的序列
				//fmt.Println("chr complete:", currentRecord.Chromosome)
				records = append(records, currentRecord)
				sequenceBuilder.Reset() // 重置 Builder
			}

			parts := strings.Fields(line)
			chromosome := parts[0]

			// 创建新的记录
			currentRecord = &FASTARecord{
				Chromosome: chromosome[1:], // 提取染色体标识符，去掉前缀 ">"
				Sequence:   "",
			}
		} else if currentRecord != nil {
			// 如果是序列数据行，将其追加到当前记录的序列中
			//currentRecord.Sequence += line //效率比较低，改成strings.Builder
			// 如果是序列数据行，将其追加到 Builder 中
			sequenceBuilder.WriteString(line)
		}
	}

	// 将最后一个记录添加到结果集中
	if currentRecord != nil {
		records = append(records, currentRecord)
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return records, nil
}

// FindCGPositions 查找序列中所有 "CG" 的出现位置
func (f *FASTARecord) FindCGPositions() []int {
	var positions []int
	sequenceLength := len(f.Sequence)

	// 遍历序列，寻找 "CG" 模式
	for i := 0; i < sequenceLength-1; i++ {
		if f.Sequence[i] == 'C' && f.Sequence[i+1] == 'G' {
			positions = append(positions, i+1) //从1开始计数
		}
	}

	return positions
}
