package simplelib

import (
	"fmt"
	"simplelib/build/go/Example/simplelib"
)

func main() {

	simpleClass := simplelib.NewSimpleClass()
	result := simpleClass.Hello()
	fmt.Println(result)
}
