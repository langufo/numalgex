#ifndef BNDRYPROP_HPP
#define BNDRYPROP_HPP

/**
 * An unsigned integer type for a bitfield to flag the sides of the
 * boundary that share some property. It can be manipulated using the static
 * members TOPBNDRY, LEFTBNDRY, RIGHTBNDRY, BOTTOMBNDRY defined in this
 * class.
 */
typedef unsigned int BndryProp;

/**
 * Bitmask to set the bit associated to the top side of the boundary.
 */
static const BndryProp TOPBNDRY = 1;

/**
 * Bitmask to set the bit associated to the top side of the boundary.
 */
static const BndryProp LEFTBNDRY = 2;

/**
 * Bitmask to set the bit associated to the right side of the boundary.
 */
static const BndryProp RIGHTBNDRY = 4;

/**
 * Bitmask to set the bit associated to the bottom side of the boundary.
 */
static const BndryProp BOTTOMBNDRY = 8;

#endif
