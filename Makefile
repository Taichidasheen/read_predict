export CGO_CFLAGS="-I/home/yangjianLab/yangwen/software/libtensorflow/include"
export CGO_LDFLAGS="-L/home/yangjianLab/yangwen/software/libtensorflow/lib -ltensorflow"


build:
	go build -o ./bin/MeLoDe_Read ./cmd/cmd.go