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
}
