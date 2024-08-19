export CGO_CFLAGS="-I/home/yangjianLab/yangwen/software/libtensorflow/include"

build:
	go build -o ./bin/MeLoDe_Read ./cmd/cmd.go