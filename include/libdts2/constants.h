#pragma once
#ifndef LIB_DTS2_CONSTANTS_H
#define LIB_DTS2_CONSTANTS_H

#define LIB_DTS2_NAMESPACE dts2
#define LIB_DTS2_ORIGIN_X int(0)
#define LIB_DTS2_ORIGIN_Y int(0)
#define LIB_DTS2_ORIGIN_Z int(0)
#define LIB_DTS2_ORIGIN Point_3(0, 0, 0)

enum class AuxPointSelector : int {
	UPPER_0=0,
	UPPER_1=1,
	UPPER_2=2,
	LOWER=3
};

#endif
