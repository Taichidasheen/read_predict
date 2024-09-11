#### 1. install tensorflow for c ####
https://tensorflow.google.cn/install/lang_c

#### 2. update the Makefile ####
Update the `Makefile` to include the TensorFlow headers
```makefile
export CGO_CFLAGS="-I/your/tensorflow/install/path/include"
```

Replace /your/tensorflow/install/path with the actual path where you installed TensorFlow.


#### 3. build executable file ####
Run the `make` command to compile and build the executable file.