extern crate bindgen;
extern crate cmake;

use cmake::Config;
use std::{env, path::PathBuf};

fn main() {
    let dst = Config::new("../../")
        .generator("Ninja")
        .define("CMAKE_BUILD_TYPE", "Release")
        .build();

    let bindings = bindgen::Builder::default()
        .header("wrapper.h")
        .clang_arg("-x")
        .clang_arg("c++")
        .clang_arg("-std=c++20")
        .clang_arg("-I../../omap")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .expect("Failed to generate bindings");
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Failed to write bindings");

    // To compile this, you need openssl static libs.
    // They exist on most distros, but you probably need to install them.
    //
    println!("cargo:rustc-link-search=native=/usr/lib");
    println!("cargo:rustc-link-search=native=/usr/lib/gcc/x86_64-linux-gnu/11/");
    println!("cargo:rustc-link-lib=static=stdc++");
    println!("cargo:rustc-link-lib=boost_system");
    println!("cargo:rustc-link-search=native=./BearSSL/build/");
    println!("cargo:rustc-link-lib=static=bearssl");
    println!("cargo:rustc-link-lib=crypto");
    println!("cargo:rustc-link-lib=static=gomp");

    println!(
        "cargo:rustc-link-search=native={}/build/omap",
        dst.display()
    );
    println!("cargo:rustc-link-lib=static=common");
    println!("cargo:rustc-link-lib=static=oraminterface");
}
