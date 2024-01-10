package worker

import "testing"

func Test_genMMTag(t *testing.T) {

	readSeqString := "AGTCTAGACTCCGTAATTACTCGCCTAG"
	aux, err := genMMTag(readSeqString)
	if err != nil {
		t.Errorf("err:%v", err)
		return
	}
	t.Logf("tag:%v, value:%v", aux.Tag(), aux.String())

}
