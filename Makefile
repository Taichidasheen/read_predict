export CGO_CFLAGS="-I/home/yangjianLab/yangwen/software/libtensorflow/include"

build:
	go build -o ./bin/read_predict ./cmd/cmd.go