#ifndef BNDRYPROP_HPP
#define BNDRYPROP_HPP

/**
 * Alias di un tipo intero per rappresentare un bitfield che
 * specifica quali lati di un bordo condividono una certa
 * proprietà
 */
typedef unsigned int BndryProp;

/**
 * Costante che imposta il bit relativo al bordo superiore
 */
const BndryProp TOPBNDRY = 1;

/**
 * Costante che imposta il bit relativo al bordo sinistro
 */
const BndryProp LEFTBNDRY = 2;

/**
 * Costante che imposta il bit relativo al bordo destro
 */
const BndryProp RIGHTBNDRY = 4;

/**
 * Costante che imposta il bit relativo al bordo inferiore
 */
const BndryProp BOTTOMBNDRY = 8;

#endif
