package record_tag

import "testing"

func Test_genMMTag(t *testing.T) {

	readSeqString := "AGTCTAGACTCCGTAATTACTCGCCTAG"
	aux, err := GenMMTag(readSeqString)
	if err != nil {
		t.Errorf("err:%v", err)
		return
	}
	t.Logf("tag:%v, value:%v", aux.Tag(), aux.String())

}
