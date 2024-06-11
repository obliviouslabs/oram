#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

// Create a wrapper around the generated bindings for omap to make it more idiomatic
// use macro to access the bindings for different key and value sizes

use std::ffi::c_void;

fn getConstCVoidPtr<T>(valptr: &T) -> *const c_void {
    let ptr = valptr as *const _ as *const c_void;
    ptr
}

fn getMutCVoidPtr<T>(valptr: &mut T) -> *mut c_void {
    let ptr = valptr as *mut _ as *mut c_void;
    ptr
}

pub enum OMapBinding {
    OMapBinding_8_8(OMapBinding_8_8),
    OMapBinding_20_32(OMapBinding_20_32),
}

use std::marker::PhantomData;
pub struct OMap<K, V> {
    inner: OMapBinding,
    _marker: PhantomData<(K, V)>,
}

// Define a macro to implement OMap for different key-value combinations
macro_rules! impl_omap {
    ($k:ty, $v:ty, $variant:ident) => {
        impl OMap<$k, $v> {
            pub fn new() -> Self {
                Self {
                    inner: unsafe { OMapBinding::$variant($variant::new()) },
                    _marker: PhantomData,
                }
            }

            pub fn init_empty(&mut self, sz: u32) {
                if let OMapBinding::$variant(omap) = &mut self.inner {
                    unsafe {
                        (*omap).InitEmpty(sz);
                    }
                }
            }

            pub fn init_empty_external(&mut self, sz: u32, cache_in_byte: u64) {
                if let OMapBinding::$variant(omap) = &mut self.inner {
                    unsafe {
                        (*omap).InitEmptyExternal(sz, cache_in_byte);
                    }
                }
            }

            pub fn insert(&mut self, key: &$k, val: &$v, isOblivious: bool) {
                if let OMapBinding::$variant(omap) = &mut self.inner {
                    unsafe {
                        if isOblivious {
                            (*omap).OInsert(getConstCVoidPtr(key), getConstCVoidPtr(val));
                        } else {
                            (*omap).Insert(getConstCVoidPtr(key), getConstCVoidPtr(val));
                        }
                    }
                }
            }

            pub fn get(&mut self, key: &$k) -> Option<$v> {
                // let mut rv: $v = 0;
                let mut rv: $v = Default::default();
                if let OMapBinding::$variant(omap) = &mut self.inner {
                    let flag =
                        unsafe { (*omap).Find(getConstCVoidPtr(key), getMutCVoidPtr(&mut rv)) };
                    if flag {
                        Some(rv)
                    } else {
                        None
                    }
                } else {
                    None
                }
            }

            pub fn erase(&mut self, key: &$k, isOblivious: bool) {
                if let OMapBinding::$variant(omap) = &mut self.inner {
                    unsafe {
                        if isOblivious {
                            (*omap).OErase(getConstCVoidPtr(key));
                        } else {
                            (*omap).Erase(getConstCVoidPtr(key));
                        }
                    }
                }
            }

            pub fn start_init(&mut self, sz: u32) {
                if let OMapBinding::$variant(omap) = &mut self.inner {
                    unsafe {
                        (*omap).StartInit(sz);
                    }
                }
            }

            pub fn finish_init(&mut self) {
                if let OMapBinding::$variant(omap) = &mut self.inner {
                    unsafe {
                        (*omap).FinishInit();
                    }
                }
            }

            pub fn start_init_external(&mut self, sz: u32, cache_in_byte: u64) {
                if let OMapBinding::$variant(omap) = &mut self.inner {
                    unsafe {
                        (*omap).StartInitExternal(sz, cache_in_byte);
                    }
                }
            }
        }
    };
}

macro_rules! drop_omap_impl {
    ($($variant:ident),*) => {
        impl<K, V> Drop for OMap<K, V> {
            fn drop(&mut self) {
                match &mut self.inner {
                    $(
                        OMapBinding::$variant(omap) => unsafe {
                            omap.Destroy();
                        },
                    )*
                }
            }
        }
    };
}

impl_omap!(u64, u64, OMapBinding_8_8);
impl_omap!([u8; 20], [u8; 32], OMapBinding_20_32);

drop_omap_impl!(OMapBinding_8_8, OMapBinding_20_32);

#[cfg(test)]
mod tests {

    use super::*;

    // rewrite the test to use the OMap struct

    #[test]
    fn test_omap_insert_and_get() {
        let sz = 100u32;
        let mut omap: OMap<u64, u64> = OMap::<u64, u64>::new();
        omap.init_empty(sz);
        omap.insert(&123u64, &456u64, false);
        omap.insert(&789u64, &101112u64, true);
        let rv1 = omap.get(&123u64);
        let rv2 = omap.get(&789u64);
        assert_eq!(rv1, Some(456u64));
        assert_eq!(rv2, Some(101112u64));
    }

    #[test]
    fn test_omap_larger_kv() {
        let sz = 100u32;
        let mut omap: OMap<[u8; 20], [u8; 32]> = OMap::<[u8; 20], [u8; 32]>::new();
        omap.init_empty(sz);
        let key1 = [1u8; 20];
        let key2 = [2u8; 20];
        let v1 = [3u8; 32];
        let v2 = [4u8; 32];
        omap.insert(&key1, &v1, false);
        omap.insert(&key2, &v2, true);
        let rv1 = omap.get(&key1);
        let rv2 = omap.get(&key2);

        assert_eq!(rv1, Some(v1));
        assert_eq!(rv2, Some(v2));
    }

    #[test]
    fn test_omap_init() {
        let sz = 1000u32;
        let mut omap: OMap<u64, u64> = OMap::<u64, u64>::new();
        omap.start_init(sz);
        omap.insert(&123u64, &456u64, false);
        omap.insert(&789u64, &101112u64, true);
        omap.finish_init();
        omap.insert(&432u64, &10u64, false);
        let rv1 = omap.get(&123u64);
        let rv2 = omap.get(&789u64);
        let rv3 = omap.get(&432u64);
        let rv4 = omap.get(&999u64);
        assert_eq!(rv1, Some(456u64));
        assert_eq!(rv2, Some(101112u64));
        assert_eq!(rv3, Some(10u64));
        assert_eq!(rv4, None);
    }

    #[test]
    fn test_omap_init_ext_mem() {
        let sz = 100000u32;
        let mut omap: OMap<u64, u64> = OMap::<u64, u64>::new();
        omap.start_init_external(sz, 2000000u64);
        omap.insert(&123u64, &456u64, false);
        omap.insert(&789u64, &101112u64, true);
        omap.finish_init();
        omap.insert(&432u64, &10u64, false);
        let rv1 = omap.get(&123u64);
        let rv2 = omap.get(&789u64);
        let rv3 = omap.get(&432u64);
        let rv4 = omap.get(&999u64);
        assert_eq!(rv1, Some(456u64));
        assert_eq!(rv2, Some(101112u64));
        assert_eq!(rv3, Some(10u64));
        assert_eq!(rv4, None);
    }

    #[test]
    fn test_omap_erase() {
        let sz = 1000u32;
        let mut omap: OMap<u64, u64> = OMap::<u64, u64>::new();
        omap.init_empty(sz);

        // insert two key-value pairs
        omap.insert(&123u64, &456u64, false);
        omap.insert(&789u64, &101112u64, true);

        // erase the first key-value pair
        omap.erase(&123u64, false);

        // check if the first key-value pair is erased
        let rv1 = omap.get(&123u64);
        assert_eq!(rv1, None);

        // check if the second key-value pair is still there
        let rv2 = omap.get(&789u64);
        assert_eq!(rv2, Some(101112u64));
    }

    #[test]
    fn test_omap_seq() {
        let sz = 100000u32;
        let mut omap: OMap<u64, u64> = OMap::<u64, u64>::new();
        omap.init_empty(sz);
        let mut keys = vec![0u64; 1000];
        let mut values = vec![0u64; 1000];
        for i in 0..1000 {
            keys[i] = i as u64;
            values[i] = i as u64 * 2;
            omap.insert(&keys[i], &values[i], false);
        }
        for i in 0..1000 {
            let rv = omap.get(&keys[i]);
            assert_eq!(rv, Some(values[i]));
        }
    }

    #[test]
    fn test_omap_random() {
        let sz = 100000u32;
        let mut omap: OMap<u64, u64> = OMap::<u64, u64>::new();
        let mut refmap: std::collections::HashMap<u64, u64> = std::collections::HashMap::new();
        omap.init_empty(sz);
        for _ in 0..sz {
            if rand::random::<bool>() {
                let key = rand::random::<u64>();
                let value = rand::random::<u64>();
                omap.insert(&key, &value, false);
                refmap.insert(key, value);
            } else if rand::random::<bool>() && !refmap.is_empty() {
                let key = *refmap
                    .keys()
                    .nth(rand::random::<usize>() % refmap.len())
                    .unwrap();
                omap.erase(&key, false);
                refmap.remove(&key);
            }
            let searchkey = rand::random::<u64>();
            let rv = omap.get(&searchkey);
            let refrv = refmap.get(&searchkey);
            if let Some(v) = refrv {
                assert_eq!(rv, Some(*v));
            } else {
                assert_eq!(rv, None);
            }
        }
    }

    // the par omaps just support u64 key and value for now and doesn't have a wrapper
    #[test]
    fn test_par_omap() {
        let sz = 10000u32;
        let mut omap: ParOMapBinding = unsafe { ParOMapBinding::new() };
        unsafe {
            omap.InitEmpty(sz, 4u32);
            let mut flags = vec![false; 4];
            let keys = vec![123u64, 456u64, 789u64, 101112u64];
            let values = vec![1u64, 2u64, 3u64, 4u64];
            omap.InsertBatch(4u32, keys.as_ptr(), values.as_ptr(), flags.as_mut_ptr());
        };
        let mut resVec = vec![0u64; 5];
        let keys = vec![789u64, 456u64, 678u64, 123u64, 456u64];
        let mut flags = vec![false; 5];
        unsafe {
            omap.FindBatch(5u32, keys.as_ptr(), resVec.as_mut_ptr(), flags.as_mut_ptr());
        };
        assert_eq!(resVec[0], 3u64);
        assert_eq!(resVec[1], 2u64);
        assert_eq!(resVec[2], 0u64);
        assert_eq!(resVec[3], 1u64);
        assert_eq!(resVec[4], 2u64);
        assert_eq!(flags[0], true);
        assert_eq!(flags[1], true);
        assert_eq!(flags[2], false);
        assert_eq!(flags[3], true);
        assert_eq!(flags[4], true);
    }

    #[test]
    fn test_par_omap_init() {
        let sz = 10000u32;
        let mut omap: ParOMapBinding = unsafe { ParOMapBinding::new() };
        unsafe {
            omap.StartInit(sz, 4u32, 4u32);
            let mut flags = vec![false; 2];
            let keys1 = vec![789u64, 123u64];
            let values1 = vec![3u64, 1u64];
            omap.InsertBatch(2u32, keys1.as_ptr(), values1.as_ptr(), flags.as_mut_ptr());
            let keys2 = vec![456u64, 101112u64];
            let values2 = vec![2u64, 4u64];
            omap.InsertBatch(2u32, keys2.as_ptr(), values2.as_ptr(), flags.as_mut_ptr());
            omap.FinishInit();
            let mut flags = vec![false; 4];
            let keys = vec![654u64, 123u64, 789u64, 678u64];
            let values = vec![1u64, 2u64, 3u64, 4u64];
            omap.InsertBatch(4u32, keys.as_ptr(), values.as_ptr(), flags.as_mut_ptr());
        };
        let mut resVec = vec![0u64; 5];
        let keys = vec![789u64, 456u64, 678u64, 123u64, 456u64];
        let mut flags = vec![false; 5];
        unsafe {
            omap.FindBatchDeferMaintain(
                5u32,
                keys.as_ptr(),
                resVec.as_mut_ptr(),
                flags.as_mut_ptr(),
            );
        };
        assert_eq!(resVec[0], 3u64);
        assert_eq!(resVec[1], 2u64);
        assert_eq!(resVec[2], 4u64);
        assert_eq!(resVec[3], 2u64);
        assert_eq!(resVec[4], 2u64);
        assert_eq!(flags[0], true);
        assert_eq!(flags[1], true);
        assert_eq!(flags[2], true);
        assert_eq!(flags[3], true);
        assert_eq!(flags[4], true);
        unsafe {
            omap.FindBatchMaintain();
        };
        let mut resVec = vec![0u64; 5];
        let keys = vec![4u64, 456u64, 678u64, 101112u64, 456u64];
        let mut flags = vec![false; 5];
        unsafe {
            omap.FindBatch(5u32, keys.as_ptr(), resVec.as_mut_ptr(), flags.as_mut_ptr());
        };
        assert_eq!(resVec[0], 0u64);
        assert_eq!(resVec[1], 2u64);
        assert_eq!(resVec[2], 4u64);
        assert_eq!(resVec[3], 4u64);
        assert_eq!(resVec[4], 2u64);
        assert_eq!(flags[0], false);
        assert_eq!(flags[1], true);
        assert_eq!(flags[2], true);
        assert_eq!(flags[3], true);
        assert_eq!(flags[4], true);
    }
}
