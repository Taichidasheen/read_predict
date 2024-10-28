package predict

import (
	"fmt"
	tf "github.com/wamuir/graft/tensorflow"
	"log"
)

func Predict(model *tf.SavedModel, inputData interface{}) ([]float32, error) {
	//log.Printf("model:%+v", model)
	//log.Printf("inputData:%v", inputData)
	/*start := time.Now()
	defer func() {
		log.Println("predict cost:", time.Since(start))
	}()*/

	tensor, err := tf.NewTensor(inputData)
	if err != nil {
		log.Printf("tf.NewTensor err: %v", err)
		return nil, err
	}

	result, err := model.Session.Run(
		map[tf.Output]*tf.Tensor{
			// python版tensorflow/keras中定义的输入层input_layer
			model.Graph.Operation("serving_default_input_1").Output(0): tensor,
		},
		[]tf.Output{
			// python版tensorflow/keras中定义的输出层output_layer
			model.Graph.Operation("StatefulPartitionedCall").Output(0),
		},
		nil,
	)
	if err != nil {
		log.Printf("model.Session.Run err:%v", err)
		return nil, err
	}
	//log.Printf("result[0].Value:%+v", result[0].Value())
	if len(result) < 1 {
		log.Printf("predict get empty result:%+v, input:%+v", result, inputData)
		return nil, fmt.Errorf("predict get empty result")
	}
	scores, ok := result[0].Value().([][]float32)
	if !ok {
		log.Printf("not expected format, result:%+v, input:%v", result, inputData)
		return nil, fmt.Errorf("not expected format")
	}
	probes := make([]float32, len(scores))
	for i, arr := range scores {
		probes[i] = arr[0]
	}
	return probes, nil
}

func BatchPredict(model *tf.SavedModel, inputs []interface{}) (interface{}, error) {
	return nil, nil
}
