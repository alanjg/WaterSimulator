#include "stdafx.h"
#include "Tests.h"
#include "LevelSet2D.h"
#include "StableFluid2D.h"

bool TestLevelSet2D()
{
	LevelSet2D lset(100, 80);
	StableFluid2D fluid(100, 80, 100, 80);
	int index = fluid.makeIndex(4, 7);
	int x, y;
	fluid.splitIndex(index, x, y);
	if (x != 4)
	{
		return false;
	}
	if (y != 7)
	{
		return false;
	}

	return true;
}

bool RunTests()
{
	return TestLevelSet2D();
}