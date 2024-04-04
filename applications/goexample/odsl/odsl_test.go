package odsl

import (
	"testing"
)

func TestORAM(t *testing.T) {
	var oram ORAMBindingSingleton
	oram.InitORAM(SwigcptrUint32_t(1000))

}
