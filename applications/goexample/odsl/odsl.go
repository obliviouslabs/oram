/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * This file is not intended to be easily readable and contains a number of
 * coding conventions designed to improve portability and efficiency. Do not make
 * changes to this file unless you know what you are doing--modify the SWIG
 * interface file instead.
 * ----------------------------------------------------------------------------- */

// source: /home/tyg/omap/oram/applications/goexample/odsl/odsl.i

package main

/*
#define intgo swig_intgo
typedef void *swig_voidp;
#cgo LDFLAGS: -lodsl -lbearssl -lstdc++
#include <stdint.h>


typedef long long intgo;
typedef unsigned long long uintgo;



typedef struct { char *p; intgo n; } _gostring_;
typedef struct { void* array; intgo len; intgo cap; } _goslice_;


typedef long long swig_type_1;
typedef long long swig_type_2;
typedef long long swig_type_3;
typedef long long swig_type_4;
typedef long long swig_type_5;
typedef long long swig_type_6;
typedef long long swig_type_7;
typedef long long swig_type_8;
typedef long long swig_type_9;
typedef long long swig_type_10;
typedef long long swig_type_11;
typedef long long swig_type_12;
typedef long long swig_type_13;
typedef long long swig_type_14;
typedef long long swig_type_15;
extern void _wrap_Swig_free_odsl_4d6b5a1bb13c6d95(uintptr_t arg1);
extern uintptr_t _wrap_Swig_malloc_odsl_4d6b5a1bb13c6d95(swig_intgo arg1);
extern void _wrap_ResetBackend_odsl_4d6b5a1bb13c6d95(swig_type_1 arg1);
extern void _wrap_DeleteBackend_odsl_4d6b5a1bb13c6d95(void);
extern void _wrap_HelloWorld_odsl_4d6b5a1bb13c6d95(swig_intgo arg1);
extern void _wrap_OMapBindingSingleton_omap_set_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, uintptr_t arg2);
extern uintptr_t _wrap_OMapBindingSingleton_omap_get_odsl_4d6b5a1bb13c6d95(uintptr_t arg1);
extern void _wrap_OMapBindingSingleton_initializer_set_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, uintptr_t arg2);
extern uintptr_t _wrap_OMapBindingSingleton_initializer_get_odsl_4d6b5a1bb13c6d95(uintptr_t arg1);
extern uintptr_t _wrap_new_OMapBindingSingleton_odsl_4d6b5a1bb13c6d95(void);
extern void _wrap_OMapBindingSingleton_InitEmpty_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_intgo arg2);
extern void _wrap_OMapBindingSingleton_InitEmptyExternal_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_intgo arg2, swig_type_2 arg3);
extern _Bool _wrap_OMapBindingSingleton_Insert_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_type_3 arg2, swig_type_4 arg3);
extern _Bool _wrap_OMapBindingSingleton_OInsert_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_type_5 arg2, swig_type_6 arg3);
extern _Bool _wrap_OMapBindingSingleton_Find_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_type_7 arg2, swig_voidp arg3);
extern _Bool _wrap_OMapBindingSingleton_Erase_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_type_8 arg2);
extern _Bool _wrap_OMapBindingSingleton_OErase_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_type_9 arg2);
extern void _wrap_OMapBindingSingleton_StartInit_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_intgo arg2);
extern void _wrap_OMapBindingSingleton_StartInitExternal_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_intgo arg2, swig_type_10 arg3);
extern void _wrap_OMapBindingSingleton_FinishInit_odsl_4d6b5a1bb13c6d95(uintptr_t arg1);
extern void _wrap_delete_OMapBindingSingleton_odsl_4d6b5a1bb13c6d95(uintptr_t arg1);
extern void _wrap_ParOMapBindingSingleton_omap_set_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, uintptr_t arg2);
extern uintptr_t _wrap_ParOMapBindingSingleton_omap_get_odsl_4d6b5a1bb13c6d95(uintptr_t arg1);
extern void _wrap_ParOMapBindingSingleton_initializer_set_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, uintptr_t arg2);
extern uintptr_t _wrap_ParOMapBindingSingleton_initializer_get_odsl_4d6b5a1bb13c6d95(uintptr_t arg1);
extern uintptr_t _wrap_new_ParOMapBindingSingleton_odsl_4d6b5a1bb13c6d95(void);
extern void _wrap_ParOMapBindingSingleton_InitEmpty_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_intgo arg2, swig_intgo arg3);
extern void _wrap_ParOMapBindingSingleton_InitEmptyExternal_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_intgo arg2, swig_intgo arg3, swig_type_11 arg4);
extern void _wrap_ParOMapBindingSingleton_InsertBatch_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_intgo arg2, swig_voidp arg3, swig_voidp arg4, swig_voidp arg5);
extern void _wrap_ParOMapBindingSingleton_FindBatch_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_intgo arg2, swig_voidp arg3, swig_voidp arg4, swig_voidp arg5);
extern void _wrap_ParOMapBindingSingleton_FindBatchDeferMaintain_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_intgo arg2, swig_voidp arg3, swig_voidp arg4, swig_voidp arg5);
extern void _wrap_ParOMapBindingSingleton_FindBatchMaintain_odsl_4d6b5a1bb13c6d95(uintptr_t arg1);
extern void _wrap_ParOMapBindingSingleton_EraseBatch_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_intgo arg2, swig_voidp arg3, swig_voidp arg4);
extern void _wrap_ParOMapBindingSingleton_StartInit_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_intgo arg2, swig_intgo arg3, swig_intgo arg4);
extern void _wrap_ParOMapBindingSingleton_StartInitExternal_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_intgo arg2, swig_intgo arg3, swig_intgo arg4, swig_type_12 arg5);
extern void _wrap_ParOMapBindingSingleton_FinishInit_odsl_4d6b5a1bb13c6d95(uintptr_t arg1);
extern void _wrap_delete_ParOMapBindingSingleton_odsl_4d6b5a1bb13c6d95(uintptr_t arg1);
extern void _wrap_ORAMBindingSingleton_oram_set_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, uintptr_t arg2);
extern uintptr_t _wrap_ORAMBindingSingleton_oram_get_odsl_4d6b5a1bb13c6d95(uintptr_t arg1);
extern uintptr_t _wrap_new_ORAMBindingSingleton_odsl_4d6b5a1bb13c6d95(void);
extern void _wrap_ORAMBindingSingleton_InitORAM_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_intgo arg2);
extern void _wrap_ORAMBindingSingleton_InitORAMExternal_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_intgo arg2, swig_type_13 arg3);
extern void _wrap_ORAMBindingSingleton_Write_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_intgo arg2, swig_type_14 arg3);
extern swig_type_15 _wrap_ORAMBindingSingleton_Read_odsl_4d6b5a1bb13c6d95(uintptr_t arg1, swig_intgo arg2);
extern void _wrap_delete_ORAMBindingSingleton_odsl_4d6b5a1bb13c6d95(uintptr_t arg1);
#undef intgo
*/
import "C"

import (
	_ "runtime/cgo"
	"sync"
	"unsafe"
)

type _ unsafe.Pointer

var Swig_escape_always_false bool
var Swig_escape_val interface{}

type _swig_fnptr *byte
type _swig_memberptr *byte

type _ sync.Mutex

func Swig_free(arg1 uintptr) {
	_swig_i_0 := arg1
	C._wrap_Swig_free_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0))
}

func Swig_malloc(arg1 int) (_swig_ret uintptr) {
	var swig_r uintptr
	_swig_i_0 := arg1
	swig_r = (uintptr)(C._wrap_Swig_malloc_odsl_4d6b5a1bb13c6d95(C.swig_intgo(_swig_i_0)))
	return swig_r
}

func ResetBackend(arg1 uint64) {
	_swig_i_0 := arg1
	C._wrap_ResetBackend_odsl_4d6b5a1bb13c6d95(C.swig_type_1(_swig_i_0))
}

func DeleteBackend() {
	C._wrap_DeleteBackend_odsl_4d6b5a1bb13c6d95()
}

func HelloWorld(arg1 uint) {
	_swig_i_0 := arg1
	C._wrap_HelloWorld_odsl_4d6b5a1bb13c6d95(C.swig_intgo(_swig_i_0))
}

type SwigcptrOMapBindingSingleton uintptr

func (p SwigcptrOMapBindingSingleton) Swigcptr() uintptr {
	return (uintptr)(p)
}

func (p SwigcptrOMapBindingSingleton) SwigIsOMapBindingSingleton() {
}

func (arg1 SwigcptrOMapBindingSingleton) SetOmap(arg2 uintptr) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	C._wrap_OMapBindingSingleton_omap_set_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.uintptr_t(_swig_i_1))
}

func (arg1 SwigcptrOMapBindingSingleton) GetOmap() (_swig_ret uintptr) {
	var swig_r uintptr
	_swig_i_0 := arg1
	swig_r = (uintptr)(C._wrap_OMapBindingSingleton_omap_get_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0)))
	return swig_r
}

func (arg1 SwigcptrOMapBindingSingleton) SetInitializer(arg2 uintptr) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	C._wrap_OMapBindingSingleton_initializer_set_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.uintptr_t(_swig_i_1))
}

func (arg1 SwigcptrOMapBindingSingleton) GetInitializer() (_swig_ret uintptr) {
	var swig_r uintptr
	_swig_i_0 := arg1
	swig_r = (uintptr)(C._wrap_OMapBindingSingleton_initializer_get_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0)))
	return swig_r
}

func NewOMapBindingSingleton() (_swig_ret OMapBindingSingleton) {
	var swig_r OMapBindingSingleton
	swig_r = (OMapBindingSingleton)(SwigcptrOMapBindingSingleton(C._wrap_new_OMapBindingSingleton_odsl_4d6b5a1bb13c6d95()))
	return swig_r
}

func (arg1 SwigcptrOMapBindingSingleton) InitEmpty(arg2 uint) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	C._wrap_OMapBindingSingleton_InitEmpty_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_intgo(_swig_i_1))
}

func (arg1 SwigcptrOMapBindingSingleton) InitEmptyExternal(arg2 uint, arg3 uint64) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	_swig_i_2 := arg3
	C._wrap_OMapBindingSingleton_InitEmptyExternal_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_intgo(_swig_i_1), C.swig_type_2(_swig_i_2))
}

func (arg1 SwigcptrOMapBindingSingleton) Insert(arg2 uint64, arg3 uint64) (_swig_ret bool) {
	var swig_r bool
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	_swig_i_2 := arg3
	swig_r = (bool)(C._wrap_OMapBindingSingleton_Insert_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_type_3(_swig_i_1), C.swig_type_4(_swig_i_2)))
	return swig_r
}

func (arg1 SwigcptrOMapBindingSingleton) OInsert(arg2 uint64, arg3 uint64) (_swig_ret bool) {
	var swig_r bool
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	_swig_i_2 := arg3
	swig_r = (bool)(C._wrap_OMapBindingSingleton_OInsert_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_type_5(_swig_i_1), C.swig_type_6(_swig_i_2)))
	return swig_r
}

func (arg1 SwigcptrOMapBindingSingleton) Find(arg2 uint64, arg3 *uint64) (_swig_ret bool) {
	var swig_r bool
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	_swig_i_2 := arg3
	swig_r = (bool)(C._wrap_OMapBindingSingleton_Find_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_type_7(_swig_i_1), C.swig_voidp(_swig_i_2)))
	return swig_r
}

func (arg1 SwigcptrOMapBindingSingleton) Erase(arg2 uint64) (_swig_ret bool) {
	var swig_r bool
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	swig_r = (bool)(C._wrap_OMapBindingSingleton_Erase_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_type_8(_swig_i_1)))
	return swig_r
}

func (arg1 SwigcptrOMapBindingSingleton) OErase(arg2 uint64) (_swig_ret bool) {
	var swig_r bool
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	swig_r = (bool)(C._wrap_OMapBindingSingleton_OErase_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_type_9(_swig_i_1)))
	return swig_r
}

func (arg1 SwigcptrOMapBindingSingleton) StartInit(arg2 uint) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	C._wrap_OMapBindingSingleton_StartInit_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_intgo(_swig_i_1))
}

func (arg1 SwigcptrOMapBindingSingleton) StartInitExternal(arg2 uint, arg3 uint64) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	_swig_i_2 := arg3
	C._wrap_OMapBindingSingleton_StartInitExternal_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_intgo(_swig_i_1), C.swig_type_10(_swig_i_2))
}

func (arg1 SwigcptrOMapBindingSingleton) FinishInit() {
	_swig_i_0 := arg1
	C._wrap_OMapBindingSingleton_FinishInit_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0))
}

func DeleteOMapBindingSingleton(arg1 OMapBindingSingleton) {
	_swig_i_0 := arg1.Swigcptr()
	C._wrap_delete_OMapBindingSingleton_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0))
}

type OMapBindingSingleton interface {
	Swigcptr() uintptr
	SwigIsOMapBindingSingleton()
	SetOmap(arg2 uintptr)
	GetOmap() (_swig_ret uintptr)
	SetInitializer(arg2 uintptr)
	GetInitializer() (_swig_ret uintptr)
	InitEmpty(arg2 uint)
	InitEmptyExternal(arg2 uint, arg3 uint64)
	Insert(arg2 uint64, arg3 uint64) (_swig_ret bool)
	OInsert(arg2 uint64, arg3 uint64) (_swig_ret bool)
	Find(arg2 uint64, arg3 *uint64) (_swig_ret bool)
	Erase(arg2 uint64) (_swig_ret bool)
	OErase(arg2 uint64) (_swig_ret bool)
	StartInit(arg2 uint)
	StartInitExternal(arg2 uint, arg3 uint64)
	FinishInit()
}

type SwigcptrParOMapBindingSingleton uintptr

func (p SwigcptrParOMapBindingSingleton) Swigcptr() uintptr {
	return (uintptr)(p)
}

func (p SwigcptrParOMapBindingSingleton) SwigIsParOMapBindingSingleton() {
}

func (arg1 SwigcptrParOMapBindingSingleton) SetOmap(arg2 uintptr) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	C._wrap_ParOMapBindingSingleton_omap_set_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.uintptr_t(_swig_i_1))
}

func (arg1 SwigcptrParOMapBindingSingleton) GetOmap() (_swig_ret uintptr) {
	var swig_r uintptr
	_swig_i_0 := arg1
	swig_r = (uintptr)(C._wrap_ParOMapBindingSingleton_omap_get_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0)))
	return swig_r
}

func (arg1 SwigcptrParOMapBindingSingleton) SetInitializer(arg2 uintptr) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	C._wrap_ParOMapBindingSingleton_initializer_set_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.uintptr_t(_swig_i_1))
}

func (arg1 SwigcptrParOMapBindingSingleton) GetInitializer() (_swig_ret uintptr) {
	var swig_r uintptr
	_swig_i_0 := arg1
	swig_r = (uintptr)(C._wrap_ParOMapBindingSingleton_initializer_get_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0)))
	return swig_r
}

func NewParOMapBindingSingleton() (_swig_ret ParOMapBindingSingleton) {
	var swig_r ParOMapBindingSingleton
	swig_r = (ParOMapBindingSingleton)(SwigcptrParOMapBindingSingleton(C._wrap_new_ParOMapBindingSingleton_odsl_4d6b5a1bb13c6d95()))
	return swig_r
}

func (arg1 SwigcptrParOMapBindingSingleton) InitEmpty(arg2 uint, arg3 uint) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	_swig_i_2 := arg3
	C._wrap_ParOMapBindingSingleton_InitEmpty_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_intgo(_swig_i_1), C.swig_intgo(_swig_i_2))
}

func (arg1 SwigcptrParOMapBindingSingleton) InitEmptyExternal(arg2 uint, arg3 uint, arg4 uint64) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	_swig_i_2 := arg3
	_swig_i_3 := arg4
	C._wrap_ParOMapBindingSingleton_InitEmptyExternal_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_intgo(_swig_i_1), C.swig_intgo(_swig_i_2), C.swig_type_11(_swig_i_3))
}

func (arg1 SwigcptrParOMapBindingSingleton) InsertBatch(arg2 uint, arg3 *uint64, arg4 *uint64, arg5 *bool) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	_swig_i_2 := arg3
	_swig_i_3 := arg4
	_swig_i_4 := arg5
	C._wrap_ParOMapBindingSingleton_InsertBatch_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_intgo(_swig_i_1), C.swig_voidp(_swig_i_2), C.swig_voidp(_swig_i_3), C.swig_voidp(_swig_i_4))
}

func (arg1 SwigcptrParOMapBindingSingleton) FindBatch(arg2 uint, arg3 *uint64, arg4 *uint64, arg5 *bool) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	_swig_i_2 := arg3
	_swig_i_3 := arg4
	_swig_i_4 := arg5
	C._wrap_ParOMapBindingSingleton_FindBatch_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_intgo(_swig_i_1), C.swig_voidp(_swig_i_2), C.swig_voidp(_swig_i_3), C.swig_voidp(_swig_i_4))
}

func (arg1 SwigcptrParOMapBindingSingleton) FindBatchDeferMaintain(arg2 uint, arg3 *uint64, arg4 *uint64, arg5 *bool) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	_swig_i_2 := arg3
	_swig_i_3 := arg4
	_swig_i_4 := arg5
	C._wrap_ParOMapBindingSingleton_FindBatchDeferMaintain_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_intgo(_swig_i_1), C.swig_voidp(_swig_i_2), C.swig_voidp(_swig_i_3), C.swig_voidp(_swig_i_4))
}

func (arg1 SwigcptrParOMapBindingSingleton) FindBatchMaintain() {
	_swig_i_0 := arg1
	C._wrap_ParOMapBindingSingleton_FindBatchMaintain_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0))
}

func (arg1 SwigcptrParOMapBindingSingleton) EraseBatch(arg2 uint, arg3 *uint64, arg4 *bool) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	_swig_i_2 := arg3
	_swig_i_3 := arg4
	C._wrap_ParOMapBindingSingleton_EraseBatch_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_intgo(_swig_i_1), C.swig_voidp(_swig_i_2), C.swig_voidp(_swig_i_3))
}

func (arg1 SwigcptrParOMapBindingSingleton) StartInit(arg2 uint, arg3 uint, arg4 uint) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	_swig_i_2 := arg3
	_swig_i_3 := arg4
	C._wrap_ParOMapBindingSingleton_StartInit_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_intgo(_swig_i_1), C.swig_intgo(_swig_i_2), C.swig_intgo(_swig_i_3))
}

func (arg1 SwigcptrParOMapBindingSingleton) StartInitExternal(arg2 uint, arg3 uint, arg4 uint, arg5 uint64) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	_swig_i_2 := arg3
	_swig_i_3 := arg4
	_swig_i_4 := arg5
	C._wrap_ParOMapBindingSingleton_StartInitExternal_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_intgo(_swig_i_1), C.swig_intgo(_swig_i_2), C.swig_intgo(_swig_i_3), C.swig_type_12(_swig_i_4))
}

func (arg1 SwigcptrParOMapBindingSingleton) FinishInit() {
	_swig_i_0 := arg1
	C._wrap_ParOMapBindingSingleton_FinishInit_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0))
}

func DeleteParOMapBindingSingleton(arg1 ParOMapBindingSingleton) {
	_swig_i_0 := arg1.Swigcptr()
	C._wrap_delete_ParOMapBindingSingleton_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0))
}

type ParOMapBindingSingleton interface {
	Swigcptr() uintptr
	SwigIsParOMapBindingSingleton()
	SetOmap(arg2 uintptr)
	GetOmap() (_swig_ret uintptr)
	SetInitializer(arg2 uintptr)
	GetInitializer() (_swig_ret uintptr)
	InitEmpty(arg2 uint, arg3 uint)
	InitEmptyExternal(arg2 uint, arg3 uint, arg4 uint64)
	InsertBatch(arg2 uint, arg3 *uint64, arg4 *uint64, arg5 *bool)
	FindBatch(arg2 uint, arg3 *uint64, arg4 *uint64, arg5 *bool)
	FindBatchDeferMaintain(arg2 uint, arg3 *uint64, arg4 *uint64, arg5 *bool)
	FindBatchMaintain()
	EraseBatch(arg2 uint, arg3 *uint64, arg4 *bool)
	StartInit(arg2 uint, arg3 uint, arg4 uint)
	StartInitExternal(arg2 uint, arg3 uint, arg4 uint, arg5 uint64)
	FinishInit()
}

type SwigcptrORAMBindingSingleton uintptr

func (p SwigcptrORAMBindingSingleton) Swigcptr() uintptr {
	return (uintptr)(p)
}

func (p SwigcptrORAMBindingSingleton) SwigIsORAMBindingSingleton() {
}

func (arg1 SwigcptrORAMBindingSingleton) SetOram(arg2 uintptr) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	C._wrap_ORAMBindingSingleton_oram_set_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.uintptr_t(_swig_i_1))
}

func (arg1 SwigcptrORAMBindingSingleton) GetOram() (_swig_ret uintptr) {
	var swig_r uintptr
	_swig_i_0 := arg1
	swig_r = (uintptr)(C._wrap_ORAMBindingSingleton_oram_get_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0)))
	return swig_r
}

func NewORAMBindingSingleton() (_swig_ret ORAMBindingSingleton) {
	var swig_r ORAMBindingSingleton
	swig_r = (ORAMBindingSingleton)(SwigcptrORAMBindingSingleton(C._wrap_new_ORAMBindingSingleton_odsl_4d6b5a1bb13c6d95()))
	return swig_r
}

func (arg1 SwigcptrORAMBindingSingleton) InitORAM(arg2 uint) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	C._wrap_ORAMBindingSingleton_InitORAM_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_intgo(_swig_i_1))
}

func (arg1 SwigcptrORAMBindingSingleton) InitORAMExternal(arg2 uint, arg3 uint64) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	_swig_i_2 := arg3
	C._wrap_ORAMBindingSingleton_InitORAMExternal_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_intgo(_swig_i_1), C.swig_type_13(_swig_i_2))
}

func (arg1 SwigcptrORAMBindingSingleton) Write(arg2 uint, arg3 uint64) {
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	_swig_i_2 := arg3
	C._wrap_ORAMBindingSingleton_Write_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_intgo(_swig_i_1), C.swig_type_14(_swig_i_2))
}

func (arg1 SwigcptrORAMBindingSingleton) Read(arg2 uint) (_swig_ret uint64) {
	var swig_r uint64
	_swig_i_0 := arg1
	_swig_i_1 := arg2
	swig_r = (uint64)(C._wrap_ORAMBindingSingleton_Read_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0), C.swig_intgo(_swig_i_1)))
	return swig_r
}

func DeleteORAMBindingSingleton(arg1 ORAMBindingSingleton) {
	_swig_i_0 := arg1.Swigcptr()
	C._wrap_delete_ORAMBindingSingleton_odsl_4d6b5a1bb13c6d95(C.uintptr_t(_swig_i_0))
}

type ORAMBindingSingleton interface {
	Swigcptr() uintptr
	SwigIsORAMBindingSingleton()
	SetOram(arg2 uintptr)
	GetOram() (_swig_ret uintptr)
	InitORAM(arg2 uint)
	InitORAMExternal(arg2 uint, arg3 uint64)
	Write(arg2 uint, arg3 uint64)
	Read(arg2 uint) (_swig_ret uint64)
}

func main() {
	var oram ORAMBindingSingleton
	oram = NewORAMBindingSingleton()
	oram.InitORAM(1000)
	oram.Write(1, 2)
	val := oram.Read(1)
	println(val)
}
