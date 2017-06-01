
import LatLonSpherical from './src/latlon-spherical'
import LatLonEllipsoidal from './src/latlon-ellipsoidal'

// merge vincenty methods into LatLonEllipsoidal
import V from './src/latlon-vincenty'
for (const prop in V) {
  LatLonEllipsoidal[prop] = V[prop]
}

import LatLonVectors from './src/latlon-vectors'
import Vector3d from './src/vector3d'
import Utm from './src/utm'
import Mgrs from './src/mgrs'
import OsGridRef from './src/osgridref'
import Dms from './src/dms'

export {
  LatLonSpherical,
  LatLonEllipsoidal,
  LatLonVectors,
  Vector3d,
  Utm,
  Mgrs,
  OsGridRef,
  Dms,
}
