#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_oram() {
        let sz = 100000u32;
        let mut oraminterface: ORAMBindingSingleton = unsafe { ORAMBindingSingleton::new() };
        unsafe {
            oraminterface.InitORAM(sz);
            oraminterface.Write(0u32, 1u64);
            oraminterface.Write(1u32, 42u64);
        };
        let rv1: u64 = unsafe { oraminterface.Read(0u32) };
        let rv2: u64 = unsafe { oraminterface.Read(1u32) };
        assert_eq!(rv1, 1);
        assert_eq!(rv2, 42);
    }

    #[test]
    fn test_omap() {
        let sz = 10000u32;
        let mut oraminterface: OMapBindingSingleton = unsafe { OMapBindingSingleton::new() };
        unsafe {
            oraminterface.InitEmpty(sz);
            oraminterface.Insert(123u64, 456u64);
            oraminterface.OInsert(789u64, 101112u64);
        };
        let mut rv1: u64 = 0;
        let mut rv2: u64 = 0;
        let flag1 = unsafe { oraminterface.Find(123u64, &mut rv1) };
        let flag2 = unsafe { oraminterface.Find(789u64, &mut rv2) };
        assert_eq!(rv1, 456u64);
        assert_eq!(rv2, 101112u64);
        assert_eq!(flag1, true);
        assert_eq!(flag2, true);
    }

    #[test]
    fn test_omap_init() {
        let sz = 100000u32;
        let mut oraminterface: OMapBindingSingleton = unsafe { OMapBindingSingleton::new() };
        unsafe {
            oraminterface.StartInit(sz);
            oraminterface.Insert(123u64, 456u64);
            oraminterface.Insert(789u64, 101112u64);
            oraminterface.FinishInit();
            oraminterface.Insert(432u64, 10u64);
        };
        let mut rv1: u64 = 0;
        let mut rv2: u64 = 0;
        let mut rv3: u64 = 0;
        let mut rv4: u64 = 0;
        let flag1 = unsafe { oraminterface.Find(123u64, &mut rv1) };
        let flag2 = unsafe { oraminterface.Find(789u64, &mut rv2) };
        let flag3 = unsafe { oraminterface.Find(432u64, &mut rv3) };
        let flag4 = unsafe { oraminterface.Find(999u64, &mut rv4) };
        assert_eq!(rv1, 456u64);
        assert_eq!(rv2, 101112u64);
        assert_eq!(rv3, 10u64);
        assert_eq!(rv4, 0u64);
        assert_eq!(flag1, true);
        assert_eq!(flag2, true);
        assert_eq!(flag3, true);
        assert_eq!(flag4, false);
    }

    #[test]
    fn test_omap_init_ext_mem() {
        let sz = 100000u32;
        let mut oraminterface: OMapBindingSingleton = unsafe { OMapBindingSingleton::new() };
        unsafe {
            oraminterface.StartInitExternal(sz, 2000000u64);
            oraminterface.Insert(123u64, 456u64);
            oraminterface.Insert(789u64, 101112u64);
            oraminterface.FinishInit();
            oraminterface.Insert(432u64, 10u64);
        };
        let mut rv1: u64 = 0;
        let mut rv2: u64 = 0;
        let mut rv3: u64 = 0;
        let mut rv4: u64 = 0;
        let flag1 = unsafe { oraminterface.Find(123u64, &mut rv1) };
        let flag2 = unsafe { oraminterface.Find(789u64, &mut rv2) };
        let flag3 = unsafe { oraminterface.Find(432u64, &mut rv3) };
        let flag4 = unsafe { oraminterface.Find(999u64, &mut rv4) };
        assert_eq!(rv1, 456u64);
        assert_eq!(rv2, 101112u64);
        assert_eq!(rv3, 10u64);
        assert_eq!(rv4, 0u64);
        assert_eq!(flag1, true);
        assert_eq!(flag2, true);
        assert_eq!(flag3, true);
        assert_eq!(flag4, false);
    }

    #[test]
    fn test_par_omap() {
        let sz = 10000u32;
        let mut oraminterface: ParOMapBindingSingleton = unsafe { ParOMapBindingSingleton::new() };
        unsafe {
            oraminterface.InitEmpty(sz, 4u32);
            let mut flags = vec![false; 4];
            let keys = vec![123u64, 456u64, 789u64, 101112u64];
            let values = vec![1u64, 2u64, 3u64, 4u64];
            oraminterface.InsertBatch(4u32, keys.as_ptr(), values.as_ptr(), flags.as_mut_ptr());
        };
        let mut resVec = vec![0u64; 5];
        let keys = vec![789u64, 456u64, 678u64, 123u64, 456u64];
        let mut flags = vec![false; 5];
        unsafe {
            oraminterface.FindBatch(5u32, keys.as_ptr(), resVec.as_mut_ptr(), flags.as_mut_ptr());
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
        let mut oraminterface: ParOMapBindingSingleton = unsafe { ParOMapBindingSingleton::new() };
        unsafe {
            oraminterface.StartInit(sz, 4u32, 4u32);
            let mut flags = vec![false; 2];
            let keys1 = vec![789u64, 123u64];
            let values1 = vec![3u64, 1u64];
            oraminterface.InsertBatch(2u32, keys1.as_ptr(), values1.as_ptr(), flags.as_mut_ptr());
            let keys2 = vec![456u64, 101112u64];
            let values2 = vec![2u64, 4u64];
            oraminterface.InsertBatch(2u32, keys2.as_ptr(), values2.as_ptr(), flags.as_mut_ptr());
            oraminterface.FinishInit();
            let mut flags = vec![false; 4];
            let keys = vec![654u64, 123u64, 789u64, 678u64];
            let values = vec![1u64, 2u64, 3u64, 4u64];
            oraminterface.InsertBatch(4u32, keys.as_ptr(), values.as_ptr(), flags.as_mut_ptr());
        };
        let mut resVec = vec![0u64; 5];
        let keys = vec![789u64, 456u64, 678u64, 123u64, 456u64];
        let mut flags = vec![false; 5];
        unsafe {
            oraminterface.FindBatchDeferMaintain(
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
            oraminterface.FindBatchMaintain();
        };
        let mut resVec = vec![0u64; 5];
        let keys = vec![4u64, 456u64, 678u64, 101112u64, 456u64];
        let mut flags = vec![false; 5];
        unsafe {
            oraminterface.FindBatch(5u32, keys.as_ptr(), resVec.as_mut_ptr(), flags.as_mut_ptr());
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
