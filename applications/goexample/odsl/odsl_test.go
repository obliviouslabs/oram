package odsl

import (
	"testing"
	"unsafe"

	"github.com/stretchr/testify/assert"
)

func TestORAM(t *testing.T) {
	var oram ORAMBinding
	oram = NewORAMBinding()
	oram.InitORAM(1000)
	oram.Write(1, 2)
	val := oram.Read(1)
	assert.Equal(t, uint64(2), val)
}

func getUintPtr(num *uint64) uintptr {
	return uintptr(unsafe.Pointer(num))
}

func TestOMap(t *testing.T) {
	var omap OMapBinding
	omap = NewOMapBinding()
	omap.InitEmpty(1000)
	key := uint64(1)
	val := uint64(2)
	omap.Insert(getUintPtr(&key), getUintPtr(&val))
	var findRes uint64
	foundFlag := omap.Find(getUintPtr(&key), getUintPtr(&findRes))
	assert.Equal(t, true, foundFlag)
	assert.Equal(t, uint64(2), findRes)
}

func TestParOMap(t *testing.T) {
	var paromap ParOMapBinding
	paromap = NewParOMapBinding()
	paromap.InitEmpty(1000, 4)
	keys := []uint64{1, 2, 3, 4, 5}
	vals := []uint64{10, 20, 30, 40, 50}
	foundFlags := []bool{false, false, false, false, false}
	paromap.InsertBatch(5, &keys[0], &vals[0], &foundFlags[0])
	for i := 0; i < 5; i++ {
		assert.Equal(t, false, foundFlags[i])
	}
	searchKeys := []uint64{3, 6, 5}
	searchVals := []uint64{0, 0, 0}
	foundFlags = []bool{false, false, false}
	paromap.FindBatch(3, &searchKeys[0], &searchVals[0], &foundFlags[0])
	assert.Equal(t, true, foundFlags[0])
	assert.Equal(t, uint64(30), searchVals[0])
	assert.Equal(t, false, foundFlags[1])
	assert.Equal(t, true, foundFlags[2])
	assert.Equal(t, uint64(50), searchVals[2])
}
