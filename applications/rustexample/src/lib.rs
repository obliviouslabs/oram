#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_oram() {
        let sz = 100000u64;
        let mut oraminterface: ORAMBindingSingleton = unsafe { ORAMBindingSingleton::new() };
        unsafe { 
            oraminterface.InitORAM(sz); 
            oraminterface.Write(0u32, 1u64);
            oraminterface.Write(1u32, 42u64);
        };
        let rv1 : u64 = unsafe { oraminterface.Read(0u32) };
        let rv2 : u64 = unsafe { oraminterface.Read(1u32) };
        assert_eq!(rv1, 1);
        assert_eq!(rv2, 42);
    }
}
