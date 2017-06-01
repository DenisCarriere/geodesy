var Dms = {};
Dms.parseDMS = function(dmsStr) {
    if (typeof dmsStr == 'number' && isFinite(dmsStr)) return Number(dmsStr);
    var dms = String(dmsStr).trim().replace(/^-/, '').replace(/[NSEW]$/i, '').split(/[^0-9.,]+/);
    if (dms[dms.length-1]=='') dms.splice(dms.length-1);
    if (dms == '') return NaN;
    var deg;
    switch (dms.length) {
        case 3:
            deg = dms[0]/1 + dms[1]/60 + dms[2]/3600;
            break;
        case 2:
            deg = dms[0]/1 + dms[1]/60;
            break;
        case 1:
            deg = dms[0];
            break;
        default:
            return NaN;
    }
    if (/^-|[WS]$/i.test(dmsStr.trim())) deg = -deg;
    return Number(deg);
};
Dms.separator = '';
Dms.toDMS = function(deg, format, dp) {
    if (isNaN(deg)) return null;
    if (format === undefined) format = 'dms';
    if (dp === undefined) {
        switch (format) {
            case 'd':    case 'deg':         dp = 4; break;
            case 'dm':   case 'deg+min':     dp = 2; break;
            case 'dms':  case 'deg+min+sec': dp = 0; break;
            default:    format = 'dms'; dp = 0;
        }
    }
    deg = Math.abs(deg);
    var dms, d, m, s;
    switch (format) {
        default:
        case 'd': case 'deg':
            d = deg.toFixed(dp);
            if (d<100) d = '0' + d;
            if (d<10) d = '0' + d;
            dms = d + '°';
            break;
        case 'dm': case 'deg+min':
            d = Math.floor(deg);
            m = ((deg*60) % 60).toFixed(dp);
            if (m == 60) { m = 0; d++; }
            d = ('000'+d).slice(-3);
            if (m<10) m = '0' + m;
            dms = d + '°'+Dms.separator + m + '′';
            break;
        case 'dms': case 'deg+min+sec':
            d = Math.floor(deg);
            m = Math.floor((deg*3600)/60) % 60;
            s = (deg*3600 % 60).toFixed(dp);
            if (s == 60) { s = (0).toFixed(dp); m++; }
            if (m == 60) { m = 0; d++; }
            d = ('000'+d).slice(-3);
            m = ('00'+m).slice(-2);
            if (s<10) s = '0' + s;
            dms = d + '°'+Dms.separator + m + '′'+Dms.separator + s + '″';
            break;
    }
    return dms;
};
Dms.toLat = function(deg, format, dp) {
    var lat = Dms.toDMS(deg, format, dp);
    return lat===null ? '–' : lat.slice(1)+Dms.separator + (deg<0 ? 'S' : 'N');
};
Dms.toLon = function(deg, format, dp) {
    var lon = Dms.toDMS(deg, format, dp);
    return lon===null ? '–' : lon+Dms.separator + (deg<0 ? 'W' : 'E');
};
Dms.toBrng = function(deg, format, dp) {
    deg = (Number(deg)+360) % 360;
    var brng =  Dms.toDMS(deg, format, dp);
    return brng===null ? '–' : brng.replace('360', '0');
};
Dms.compassPoint = function(bearing, precision) {
    if (precision === undefined) precision = 3;
    bearing = ((bearing%360)+360)%360;
    var cardinals = {
        1: [ 'N', 'E', 'S', 'W' ],
        2: [ 'N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW' ],
        3: [ 'N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW' ],
    };
    var points = cardinals[precision].length;
    var cardinal = cardinals[precision][Math.round(bearing/360*points)%points];
    return cardinal;
};

function LatLon(lat, lon) {
    if (!(this instanceof LatLon)) return new LatLon(lat, lon);
    this.lat = Number(lat);
    this.lon = Number(lon);
}
LatLon.prototype.distanceTo = function(point, radius) {
    if (!(point instanceof LatLon)) throw new TypeError('point is not LatLon object');
    radius = (radius === undefined) ? 6371e3 : Number(radius);
    var R = radius;
    var φ1 = this.lat.toRadians(),  λ1 = this.lon.toRadians();
    var φ2 = point.lat.toRadians(), λ2 = point.lon.toRadians();
    var Δφ = φ2 - φ1;
    var Δλ = λ2 - λ1;
    var a = Math.sin(Δφ/2) * Math.sin(Δφ/2)
          + Math.cos(φ1) * Math.cos(φ2)
          * Math.sin(Δλ/2) * Math.sin(Δλ/2);
    var c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a));
    var d = R * c;
    return d;
};
LatLon.prototype.bearingTo = function(point) {
    if (!(point instanceof LatLon)) throw new TypeError('point is not LatLon object');
    var φ1 = this.lat.toRadians(), φ2 = point.lat.toRadians();
    var Δλ = (point.lon-this.lon).toRadians();
    var y = Math.sin(Δλ) * Math.cos(φ2);
    var x = Math.cos(φ1)*Math.sin(φ2) -
            Math.sin(φ1)*Math.cos(φ2)*Math.cos(Δλ);
    var θ = Math.atan2(y, x);
    return (θ.toDegrees()+360) % 360;
};
LatLon.prototype.finalBearingTo = function(point) {
    if (!(point instanceof LatLon)) throw new TypeError('point is not LatLon object');
    return ( point.bearingTo(this)+180 ) % 360;
};
LatLon.prototype.midpointTo = function(point) {
    if (!(point instanceof LatLon)) throw new TypeError('point is not LatLon object');
    var φ1 = this.lat.toRadians(), λ1 = this.lon.toRadians();
    var φ2 = point.lat.toRadians();
    var Δλ = (point.lon-this.lon).toRadians();
    var Bx = Math.cos(φ2) * Math.cos(Δλ);
    var By = Math.cos(φ2) * Math.sin(Δλ);
    var x = Math.sqrt((Math.cos(φ1) + Bx) * (Math.cos(φ1) + Bx) + By * By);
    var y = Math.sin(φ1) + Math.sin(φ2);
    var φ3 = Math.atan2(y, x);
    var λ3 = λ1 + Math.atan2(By, Math.cos(φ1) + Bx);
    return new LatLon(φ3.toDegrees(), (λ3.toDegrees()+540)%360-180);
};
LatLon.prototype.intermediatePointTo = function(point, fraction) {
    if (!(point instanceof LatLon)) throw new TypeError('point is not LatLon object');
    var φ1 = this.lat.toRadians(), λ1 = this.lon.toRadians();
    var φ2 = point.lat.toRadians(), λ2 = point.lon.toRadians();
    var sinφ1 = Math.sin(φ1), cosφ1 = Math.cos(φ1), sinλ1 = Math.sin(λ1), cosλ1 = Math.cos(λ1);
    var sinφ2 = Math.sin(φ2), cosφ2 = Math.cos(φ2), sinλ2 = Math.sin(λ2), cosλ2 = Math.cos(λ2);
    var Δφ = φ2 - φ1;
    var Δλ = λ2 - λ1;
    var a = Math.sin(Δφ/2) * Math.sin(Δφ/2)
        + Math.cos(φ1) * Math.cos(φ2) * Math.sin(Δλ/2) * Math.sin(Δλ/2);
    var δ = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a));
    var A = Math.sin((1-fraction)*δ) / Math.sin(δ);
    var B = Math.sin(fraction*δ) / Math.sin(δ);
    var x = A * cosφ1 * cosλ1 + B * cosφ2 * cosλ2;
    var y = A * cosφ1 * sinλ1 + B * cosφ2 * sinλ2;
    var z = A * sinφ1 + B * sinφ2;
    var φ3 = Math.atan2(z, Math.sqrt(x*x + y*y));
    var λ3 = Math.atan2(y, x);
    return new LatLon(φ3.toDegrees(), (λ3.toDegrees()+540)%360-180);
};
LatLon.prototype.destinationPoint = function(distance, bearing, radius) {
    radius = (radius === undefined) ? 6371e3 : Number(radius);
    var δ = Number(distance) / radius;
    var θ = Number(bearing).toRadians();
    var φ1 = this.lat.toRadians();
    var λ1 = this.lon.toRadians();
    var sinφ1 = Math.sin(φ1), cosφ1 = Math.cos(φ1);
    var sinδ = Math.sin(δ), cosδ = Math.cos(δ);
    var sinθ = Math.sin(θ), cosθ = Math.cos(θ);
    var sinφ2 = sinφ1*cosδ + cosφ1*sinδ*cosθ;
    var φ2 = Math.asin(sinφ2);
    var y = sinθ * sinδ * cosφ1;
    var x = cosδ - sinφ1 * sinφ2;
    var λ2 = λ1 + Math.atan2(y, x);
    return new LatLon(φ2.toDegrees(), (λ2.toDegrees()+540)%360-180);
};
LatLon.intersection = function(p1, brng1, p2, brng2) {
    if (!(p1 instanceof LatLon)) throw new TypeError('p1 is not LatLon object');
    if (!(p2 instanceof LatLon)) throw new TypeError('p2 is not LatLon object');
    var φ1 = p1.lat.toRadians(), λ1 = p1.lon.toRadians();
    var φ2 = p2.lat.toRadians(), λ2 = p2.lon.toRadians();
    var θ13 = Number(brng1).toRadians(), θ23 = Number(brng2).toRadians();
    var Δφ = φ2-φ1, Δλ = λ2-λ1;
    var δ12 = 2*Math.asin( Math.sqrt( Math.sin(Δφ/2)*Math.sin(Δφ/2)
        + Math.cos(φ1)*Math.cos(φ2)*Math.sin(Δλ/2)*Math.sin(Δλ/2) ) );
    if (δ12 == 0) return null;
    var θa = Math.acos( ( Math.sin(φ2) - Math.sin(φ1)*Math.cos(δ12) ) / ( Math.sin(δ12)*Math.cos(φ1) ) );
    if (isNaN(θa)) θa = 0;
    var θb = Math.acos( ( Math.sin(φ1) - Math.sin(φ2)*Math.cos(δ12) ) / ( Math.sin(δ12)*Math.cos(φ2) ) );
    var θ12 = Math.sin(λ2-λ1)>0 ? θa : 2*Math.PI-θa;
    var θ21 = Math.sin(λ2-λ1)>0 ? 2*Math.PI-θb : θb;
    var α1 = θ13 - θ12;
    var α2 = θ21 - θ23;
    if (Math.sin(α1)==0 && Math.sin(α2)==0) return null;
    if (Math.sin(α1)*Math.sin(α2) < 0) return null;
    var α3 = Math.acos( -Math.cos(α1)*Math.cos(α2) + Math.sin(α1)*Math.sin(α2)*Math.cos(δ12) );
    var δ13 = Math.atan2( Math.sin(δ12)*Math.sin(α1)*Math.sin(α2), Math.cos(α2)+Math.cos(α1)*Math.cos(α3) );
    var φ3 = Math.asin( Math.sin(φ1)*Math.cos(δ13) + Math.cos(φ1)*Math.sin(δ13)*Math.cos(θ13) );
    var Δλ13 = Math.atan2( Math.sin(θ13)*Math.sin(δ13)*Math.cos(φ1), Math.cos(δ13)-Math.sin(φ1)*Math.sin(φ3) );
    var λ3 = λ1 + Δλ13;
    return new LatLon(φ3.toDegrees(), (λ3.toDegrees()+540)%360-180);
};
LatLon.prototype.crossTrackDistanceTo = function(pathStart, pathEnd, radius) {
    if (!(pathStart instanceof LatLon)) throw new TypeError('pathStart is not LatLon object');
    if (!(pathEnd instanceof LatLon)) throw new TypeError('pathEnd is not LatLon object');
    var R = (radius === undefined) ? 6371e3 : Number(radius);
    var δ13 = pathStart.distanceTo(this, R) / R;
    var θ13 = pathStart.bearingTo(this).toRadians();
    var θ12 = pathStart.bearingTo(pathEnd).toRadians();
    var δxt = Math.asin(Math.sin(δ13) * Math.sin(θ13-θ12));
    return δxt * R;
};
LatLon.prototype.alongTrackDistanceTo = function(pathStart, pathEnd, radius) {
    if (!(pathStart instanceof LatLon)) throw new TypeError('pathStart is not LatLon object');
    if (!(pathEnd instanceof LatLon)) throw new TypeError('pathEnd is not LatLon object');
    var R = (radius === undefined) ? 6371e3 : Number(radius);
    var δ13 = pathStart.distanceTo(this, R) / R;
    var θ13 = pathStart.bearingTo(this).toRadians();
    var θ12 = pathStart.bearingTo(pathEnd).toRadians();
    var δxt = Math.asin(Math.sin(δ13) * Math.sin(θ13-θ12));
    var δat = Math.acos(Math.cos(δ13) / Math.abs(Math.cos(δxt)));
    return δat*Math.sign(Math.cos(θ12-θ13)) * R;
};
LatLon.prototype.maxLatitude = function(bearing) {
    var θ = Number(bearing).toRadians();
    var φ = this.lat.toRadians();
    var φMax = Math.acos(Math.abs(Math.sin(θ)*Math.cos(φ)));
    return φMax.toDegrees();
};
LatLon.crossingParallels = function(point1, point2, latitude) {
    var φ = Number(latitude).toRadians();
    var φ1 = point1.lat.toRadians();
    var λ1 = point1.lon.toRadians();
    var φ2 = point2.lat.toRadians();
    var λ2 = point2.lon.toRadians();
    var Δλ = λ2 - λ1;
    var x = Math.sin(φ1) * Math.cos(φ2) * Math.cos(φ) * Math.sin(Δλ);
    var y = Math.sin(φ1) * Math.cos(φ2) * Math.cos(φ) * Math.cos(Δλ) - Math.cos(φ1) * Math.sin(φ2) * Math.cos(φ);
    var z = Math.cos(φ1) * Math.cos(φ2) * Math.sin(φ) * Math.sin(Δλ);
    if (z*z > x*x + y*y) return null;
    var λm = Math.atan2(-y, x);
    var Δλi = Math.acos(z / Math.sqrt(x*x+y*y));
    var λi1 = λ1 + λm - Δλi;
    var λi2 = λ1 + λm + Δλi;
    return { lon1: (λi1.toDegrees()+540)%360-180, lon2: (λi2.toDegrees()+540)%360-180 };
};
LatLon.prototype.rhumbDistanceTo = function(point, radius) {
    if (!(point instanceof LatLon)) throw new TypeError('point is not LatLon object');
    radius = (radius === undefined) ? 6371e3 : Number(radius);
    var R = radius;
    var φ1 = this.lat.toRadians(), φ2 = point.lat.toRadians();
    var Δφ = φ2 - φ1;
    var Δλ = Math.abs(point.lon-this.lon).toRadians();
    if (Δλ > Math.PI) Δλ -= 2*Math.PI;
    var Δψ = Math.log(Math.tan(φ2/2+Math.PI/4)/Math.tan(φ1/2+Math.PI/4));
    var q = Math.abs(Δψ) > 10e-12 ? Δφ/Δψ : Math.cos(φ1);
    var δ = Math.sqrt(Δφ*Δφ + q*q*Δλ*Δλ);
    var dist = δ * R;
    return dist;
};
LatLon.prototype.rhumbBearingTo = function(point) {
    if (!(point instanceof LatLon)) throw new TypeError('point is not LatLon object');
    var φ1 = this.lat.toRadians(), φ2 = point.lat.toRadians();
    var Δλ = (point.lon-this.lon).toRadians();
    if (Δλ >  Math.PI) Δλ -= 2*Math.PI;
    if (Δλ < -Math.PI) Δλ += 2*Math.PI;
    var Δψ = Math.log(Math.tan(φ2/2+Math.PI/4)/Math.tan(φ1/2+Math.PI/4));
    var θ = Math.atan2(Δλ, Δψ);
    return (θ.toDegrees()+360) % 360;
};
LatLon.prototype.rhumbDestinationPoint = function(distance, bearing, radius) {
    radius = (radius === undefined) ? 6371e3 : Number(radius);
    var δ = Number(distance) / radius;
    var φ1 = this.lat.toRadians(), λ1 = this.lon.toRadians();
    var θ = Number(bearing).toRadians();
    var Δφ = δ * Math.cos(θ);
    var φ2 = φ1 + Δφ;
    if (Math.abs(φ2) > Math.PI/2) φ2 = φ2>0 ? Math.PI-φ2 : -Math.PI-φ2;
    var Δψ = Math.log(Math.tan(φ2/2+Math.PI/4)/Math.tan(φ1/2+Math.PI/4));
    var q = Math.abs(Δψ) > 10e-12 ? Δφ / Δψ : Math.cos(φ1);
    var Δλ = δ*Math.sin(θ)/q;
    var λ2 = λ1 + Δλ;
    return new LatLon(φ2.toDegrees(), (λ2.toDegrees()+540) % 360 - 180);
};
LatLon.prototype.rhumbMidpointTo = function(point) {
    if (!(point instanceof LatLon)) throw new TypeError('point is not LatLon object');
    var φ1 = this.lat.toRadians(), λ1 = this.lon.toRadians();
    var φ2 = point.lat.toRadians(), λ2 = point.lon.toRadians();
    if (Math.abs(λ2-λ1) > Math.PI) λ1 += 2*Math.PI;
    var φ3 = (φ1+φ2)/2;
    var f1 = Math.tan(Math.PI/4 + φ1/2);
    var f2 = Math.tan(Math.PI/4 + φ2/2);
    var f3 = Math.tan(Math.PI/4 + φ3/2);
    var λ3 = ( (λ2-λ1)*Math.log(f3) + λ1*Math.log(f2) - λ2*Math.log(f1) ) / Math.log(f2/f1);
    if (!isFinite(λ3)) λ3 = (λ1+λ2)/2;
    var p = LatLon(φ3.toDegrees(), (λ3.toDegrees()+540)%360-180);
    return p;
};
LatLon.areaOf = function(polygon, radius) {
    var R = (radius === undefined) ? 6371e3 : Number(radius);
    var closed = polygon[0].equals(polygon[polygon.length-1]);
    if (!closed) polygon.push(polygon[0]);
    var nVertices = polygon.length - 1;
    var S = 0;
    for (var v=0; v<nVertices; v++) {
        var φ1 = polygon[v].lat.toRadians();
        var φ2 = polygon[v+1].lat.toRadians();
        var Δλ = (polygon[v+1].lon - polygon[v].lon).toRadians();
        var E = 2 * Math.atan2(Math.tan(Δλ/2) * (Math.tan(φ1/2)+Math.tan(φ2/2)), 1 + Math.tan(φ1/2)*Math.tan(φ2/2));
        S += E;
    }
    if (isPoleEnclosedBy(polygon)) S = Math.abs(S) - 2*Math.PI;
    var A = Math.abs(S * R*R);
    if (!closed) polygon.pop();
    return A;
    function isPoleEnclosedBy(polygon) {
        var ΣΔ = 0;
        var prevBrng = polygon[0].bearingTo(polygon[1]);
        for (var v=0; v<polygon.length-1; v++) {
            var initBrng = polygon[v].bearingTo(polygon[v+1]);
            var finalBrng = polygon[v].finalBearingTo(polygon[v+1]);
            ΣΔ += (initBrng - prevBrng + 540) % 360 - 180;
            ΣΔ += (finalBrng - initBrng + 540) % 360 - 180;
            prevBrng = finalBrng;
        }
        var initBrng = polygon[0].bearingTo(polygon[1]);
        ΣΔ += (initBrng - prevBrng + 540) % 360 - 180;
        var enclosed = Math.abs(ΣΔ) < 90;
        return enclosed;
    }
};
LatLon.prototype.equals = function(point) {
    if (!(point instanceof LatLon)) throw new TypeError('point is not LatLon object');
    if (this.lat != point.lat) return false;
    if (this.lon != point.lon) return false;
    return true;
};
LatLon.prototype.toString = function(format, dp) {
    return Dms.toLat(this.lat, format, dp) + ', ' + Dms.toLon(this.lon, format, dp);
};
if (Number.prototype.toRadians === undefined) {
    Number.prototype.toRadians = function() { return this * Math.PI / 180; };
}
if (Number.prototype.toDegrees === undefined) {
    Number.prototype.toDegrees = function() { return this * 180 / Math.PI; };
}

function Vector3d(x, y, z) {
    if (!(this instanceof Vector3d)) return new Vector3d(x, y, z);
    this.x = Number(x);
    this.y = Number(y);
    this.z = Number(z);
}
Vector3d.prototype.plus = function(v) {
    if (!(v instanceof Vector3d)) throw new TypeError('v is not Vector3d object');
    return new Vector3d(this.x + v.x, this.y + v.y, this.z + v.z);
};
Vector3d.prototype.minus = function(v) {
    if (!(v instanceof Vector3d)) throw new TypeError('v is not Vector3d object');
    return new Vector3d(this.x - v.x, this.y - v.y, this.z - v.z);
};
Vector3d.prototype.times = function(x) {
    x = Number(x);
    return new Vector3d(this.x * x, this.y * x, this.z * x);
};
Vector3d.prototype.dividedBy = function(x) {
    x = Number(x);
    return new Vector3d(this.x / x, this.y / x, this.z / x);
};
Vector3d.prototype.dot = function(v) {
    if (!(v instanceof Vector3d)) throw new TypeError('v is not Vector3d object');
    return this.x*v.x + this.y*v.y + this.z*v.z;
};
Vector3d.prototype.cross = function(v) {
    if (!(v instanceof Vector3d)) throw new TypeError('v is not Vector3d object');
    var x = this.y*v.z - this.z*v.y;
    var y = this.z*v.x - this.x*v.z;
    var z = this.x*v.y - this.y*v.x;
    return new Vector3d(x, y, z);
};
Vector3d.prototype.negate = function() {
    return new Vector3d(-this.x, -this.y, -this.z);
};
Vector3d.prototype.length = function() {
    return Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z);
};
Vector3d.prototype.unit = function() {
    var norm = this.length();
    if (norm == 1) return this;
    if (norm == 0) return this;
    var x = this.x/norm;
    var y = this.y/norm;
    var z = this.z/norm;
    return new Vector3d(x, y, z);
};
Vector3d.prototype.angleTo = function(v, n) {
    if (!(v instanceof Vector3d)) throw new TypeError('v is not Vector3d object');
    if (!(n instanceof Vector3d || n == undefined)) throw new TypeError('n is not Vector3d object');
    var sign = n==undefined ? 1 : Math.sign(this.cross(v).dot(n));
    var sinθ = this.cross(v).length() * sign;
    var cosθ = this.dot(v);
    return Math.atan2(sinθ, cosθ);
};
Vector3d.prototype.rotateAround = function(axis, theta) {
    if (!(axis instanceof Vector3d)) throw new TypeError('axis is not Vector3d object');
    var p1 = this.unit();
    var p = [ p1.x, p1.y, p1.z ];
    var a = axis.unit();
    var s = Math.sin(theta);
    var c = Math.cos(theta);
    var q = [
        [ a.x*a.x*(1-c) + c,     a.x*a.y*(1-c) - a.z*s, a.x*a.z*(1-c) + a.y*s ],
        [ a.y*a.x*(1-c) + a.z*s, a.y*a.y*(1-c) + c,     a.y*a.z*(1-c) - a.x*s ],
        [ a.z*a.x*(1-c) - a.y*s, a.z*a.y*(1-c) + a.x*s, a.z*a.z*(1-c) + c     ],
    ];
    var qp = [ 0, 0, 0 ];
    for (var i=0; i<3; i++) {
        for (var j=0; j<3; j++) {
            qp[i] += q[i][j] * p[j];
        }
    }
    var p2 = new Vector3d(qp[0], qp[1], qp[2]);
    return p2;
};
Vector3d.prototype.toString = function(precision) {
    var p = (precision === undefined) ? 3 : Number(precision);
    var str = '[' + this.x.toFixed(p) + ',' + this.y.toFixed(p) + ',' + this.z.toFixed(p) + ']';
    return str;
};
if (Math.sign === undefined) {
    Math.sign = function(x) {
        x = +x;
        if (x === 0 || isNaN(x)) return x;
        return x > 0 ? 1 : -1;
    };
}

function LatLon$2(lat, lon, datum) {
    if (!(this instanceof LatLon$2)) return new LatLon$2(lat, lon, datum);
    if (datum === undefined) datum = LatLon$2.datum.WGS84;
    this.lat = Number(lat);
    this.lon = Number(lon);
    this.datum = datum;
}
LatLon$2.ellipsoid = {
    WGS84:         { a: 6378137,     b: 6356752.314245, f: 1/298.257223563 },
    Airy1830:      { a: 6377563.396, b: 6356256.909,    f: 1/299.3249646   },
    AiryModified:  { a: 6377340.189, b: 6356034.448,    f: 1/299.3249646   },
    Bessel1841:    { a: 6377397.155, b: 6356078.962818, f: 1/299.1528128   },
    Clarke1866:    { a: 6378206.4,   b: 6356583.8,      f: 1/294.978698214 },
    Clarke1880IGN: { a: 6378249.2,   b: 6356515.0,      f: 1/293.466021294 },
    GRS80:         { a: 6378137,     b: 6356752.314140, f: 1/298.257222101 },
    Intl1924:      { a: 6378388,     b: 6356911.946,    f: 1/297           },
    WGS72:         { a: 6378135,     b: 6356750.5,      f: 1/298.26        },
};
LatLon$2.datum = {
    ED50:       { ellipsoid: LatLon$2.ellipsoid.Intl1924,      transform: [   89.5,    93.8,    123.1,    -1.2,     0.0,     0.0,     0.156  ] },
    Irl1975:    { ellipsoid: LatLon$2.ellipsoid.AiryModified,  transform: [ -482.530, 130.596, -564.557,  -8.150,  -1.042,  -0.214,  -0.631  ] },
    NAD27:      { ellipsoid: LatLon$2.ellipsoid.Clarke1866,    transform: [    8,    -160,     -176,       0,       0,       0,       0      ] },
    NAD83:      { ellipsoid: LatLon$2.ellipsoid.GRS80,         transform: [    1.004,  -1.910,   -0.515,  -0.0015,  0.0267,  0.00034, 0.011  ] },
    NTF:        { ellipsoid: LatLon$2.ellipsoid.Clarke1880IGN, transform: [  168,      60,     -320,       0,       0,       0,       0      ] },
    OSGB36:     { ellipsoid: LatLon$2.ellipsoid.Airy1830,      transform: [ -446.448, 125.157, -542.060,  20.4894, -0.1502, -0.2470, -0.8421 ] },
    Potsdam:    { ellipsoid: LatLon$2.ellipsoid.Bessel1841,    transform: [ -582,    -105,     -414,      -8.3,     1.04,    0.35,   -3.08   ] },
    TokyoJapan: { ellipsoid: LatLon$2.ellipsoid.Bessel1841,    transform: [  148,    -507,     -685,       0,       0,       0,       0      ] },
    WGS72:      { ellipsoid: LatLon$2.ellipsoid.WGS72,         transform: [    0,       0,     -4.5,      -0.22,    0,       0,       0.554  ] },
    WGS84:      { ellipsoid: LatLon$2.ellipsoid.WGS84,         transform: [    0.0,     0.0,      0.0,     0.0,     0.0,     0.0,     0.0    ] },
};
LatLon$2.prototype.convertDatum = function(toDatum) {
    var oldLatLon = this;
    var transform = null;
    if (oldLatLon.datum == LatLon$2.datum.WGS84) {
        transform = toDatum.transform;
    }
    if (toDatum == LatLon$2.datum.WGS84) {
        transform = [];
        for (var p=0; p<7; p++) transform[p] = -oldLatLon.datum.transform[p];
    }
    if (transform == null) {
        oldLatLon = this.convertDatum(LatLon$2.datum.WGS84);
        transform = toDatum.transform;
    }
    var oldCartesian = oldLatLon.toCartesian();
    var newCartesian = oldCartesian.applyTransform(transform);
    var newLatLon = newCartesian.toLatLonE(toDatum);
    return newLatLon;
};
LatLon$2.prototype.toCartesian = function() {
    var φ = this.lat.toRadians(), λ = this.lon.toRadians();
    var h = 0;
    var a = this.datum.ellipsoid.a, f = this.datum.ellipsoid.f;
    var sinφ = Math.sin(φ), cosφ = Math.cos(φ);
    var sinλ = Math.sin(λ), cosλ = Math.cos(λ);
    var eSq = 2*f - f*f;
    var ν = a / Math.sqrt(1 - eSq*sinφ*sinφ);
    var x = (ν+h) * cosφ * cosλ;
    var y = (ν+h) * cosφ * sinλ;
    var z = (ν*(1-eSq)+h) * sinφ;
    var point = new Vector3d(x, y, z);
    return point;
};
Vector3d.prototype.toLatLonE = function(datum) {
    var x = this.x, y = this.y, z = this.z;
    var a = datum.ellipsoid.a, b = datum.ellipsoid.b, f = datum.ellipsoid.f;
    var e2 = 2*f - f*f;
    var ε2 = e2 / (1-e2);
    var p = Math.sqrt(x*x + y*y);
    var R = Math.sqrt(p*p + z*z);
    var tanβ = (b*z)/(a*p) * (1+ε2*b/R);
    var sinβ = tanβ / Math.sqrt(1+tanβ*tanβ);
    var cosβ = sinβ / tanβ;
    var φ = isNaN(cosβ) ? 0 : Math.atan2(z + ε2*b*sinβ*sinβ*sinβ, p - e2*a*cosβ*cosβ*cosβ);
    var λ = Math.atan2(y, x);
    var sinφ = Math.sin(φ), cosφ = Math.cos(φ);
    var ν = a/Math.sqrt(1-e2*sinφ*sinφ);
    var h = p*cosφ + z*sinφ - (a*a/ν);
    var point = new LatLon$2(φ.toDegrees(), λ.toDegrees(), datum);
    return point;
};
Vector3d.prototype.applyTransform = function(t)   {
    var x1 = this.x, y1 = this.y, z1 = this.z;
    var tx = t[0];
    var ty = t[1];
    var tz = t[2];
    var s1 = t[3]/1e6 + 1;
    var rx = (t[4]/3600).toRadians();
    var ry = (t[5]/3600).toRadians();
    var rz = (t[6]/3600).toRadians();
    var x2 = tx + x1*s1 - y1*rz + z1*ry;
    var y2 = ty + x1*rz + y1*s1 - z1*rx;
    var z2 = tz - x1*ry + y1*rx + z1*s1;
    return new Vector3d(x2, y2, z2);
};
LatLon$2.prototype.toString = function(format, dp) {
    return Dms.toLat(this.lat, format, dp) + ', ' + Dms.toLon(this.lon, format, dp);
};
if (Number.prototype.toRadians === undefined) {
    Number.prototype.toRadians = function() { return this * Math.PI / 180; };
}
if (Number.prototype.toDegrees === undefined) {
    Number.prototype.toDegrees = function() { return this * 180 / Math.PI; };
}

LatLon$2.prototype.distanceTo = function(point) {
    if (!(point instanceof LatLon$2)) throw new TypeError('point is not LatLon object');
    try {
        return Number(this.inverse(point).distance.toFixed(3));
    } catch (e) {
        return NaN;
    }
};
LatLon$2.prototype.initialBearingTo = function(point) {
    if (!(point instanceof LatLon$2)) throw new TypeError('point is not LatLon object');
    try {
        return Number(this.inverse(point).initialBearing.toFixed(9));
    } catch (e) {
        return NaN;
    }
};
LatLon$2.prototype.finalBearingTo = function(point) {
    if (!(point instanceof LatLon$2)) throw new TypeError('point is not LatLon object');
    try {
        return Number(this.inverse(point).finalBearing.toFixed(9));
    } catch (e) {
        return NaN;
    }
};
LatLon$2.prototype.destinationPoint = function(distance, initialBearing) {
    return this.direct(Number(distance), Number(initialBearing)).point;
};
LatLon$2.prototype.finalBearingOn = function(distance, initialBearing) {
    return Number(this.direct(Number(distance), Number(initialBearing)).finalBearing.toFixed(9));
};
LatLon$2.prototype.direct = function(distance, initialBearing) {
    var φ1 = this.lat.toRadians(), λ1 = this.lon.toRadians();
    var α1 = initialBearing.toRadians();
    var s = distance;
    var a = this.datum.ellipsoid.a, b = this.datum.ellipsoid.b, f = this.datum.ellipsoid.f;
    var sinα1 = Math.sin(α1);
    var cosα1 = Math.cos(α1);
    var tanU1 = (1-f) * Math.tan(φ1), cosU1 = 1 / Math.sqrt((1 + tanU1*tanU1)), sinU1 = tanU1 * cosU1;
    var σ1 = Math.atan2(tanU1, cosα1);
    var sinα = cosU1 * sinα1;
    var cosSqα = 1 - sinα*sinα;
    var uSq = cosSqα * (a*a - b*b) / (b*b);
    var A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
    var B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));
    var cos2σM, sinσ, cosσ, Δσ;
    var σ = s / (b*A), σʹ, iterations = 0;
    do {
        cos2σM = Math.cos(2*σ1 + σ);
        sinσ = Math.sin(σ);
        cosσ = Math.cos(σ);
        Δσ = B*sinσ*(cos2σM+B/4*(cosσ*(-1+2*cos2σM*cos2σM)-
            B/6*cos2σM*(-3+4*sinσ*sinσ)*(-3+4*cos2σM*cos2σM)));
        σʹ = σ;
        σ = s / (b*A) + Δσ;
    } while (Math.abs(σ-σʹ) > 1e-12 && ++iterations<100);
    if (iterations >= 100) throw new Error('Formula failed to converge');
    var x = sinU1*sinσ - cosU1*cosσ*cosα1;
    var φ2 = Math.atan2(sinU1*cosσ + cosU1*sinσ*cosα1, (1-f)*Math.sqrt(sinα*sinα + x*x));
    var λ = Math.atan2(sinσ*sinα1, cosU1*cosσ - sinU1*sinσ*cosα1);
    var C = f/16*cosSqα*(4+f*(4-3*cosSqα));
    var L = λ - (1-C) * f * sinα *
        (σ + C*sinσ*(cos2σM+C*cosσ*(-1+2*cos2σM*cos2σM)));
    var λ2 = (λ1+L+3*Math.PI)%(2*Math.PI) - Math.PI;
    var α2 = Math.atan2(sinα, -x);
    α2 = (α2 + 2*Math.PI) % (2*Math.PI);
    return {
        point:        new LatLon$2(φ2.toDegrees(), λ2.toDegrees(), this.datum),
        finalBearing: α2.toDegrees(),
        iterations:   iterations,
    };
};
LatLon$2.prototype.inverse = function(point) {
    var p1 = this, p2 = point;
    if (p1.lon == -180) p1.lon = 180;
    var φ1 = p1.lat.toRadians(), λ1 = p1.lon.toRadians();
    var φ2 = p2.lat.toRadians(), λ2 = p2.lon.toRadians();
    var a = this.datum.ellipsoid.a, b = this.datum.ellipsoid.b, f = this.datum.ellipsoid.f;
    var L = λ2 - λ1;
    var tanU1 = (1-f) * Math.tan(φ1), cosU1 = 1 / Math.sqrt((1 + tanU1*tanU1)), sinU1 = tanU1 * cosU1;
    var tanU2 = (1-f) * Math.tan(φ2), cosU2 = 1 / Math.sqrt((1 + tanU2*tanU2)), sinU2 = tanU2 * cosU2;
    var sinλ, cosλ, sinSqσ, sinσ=0, cosσ=0, σ=0, sinα, cosSqα=0, cos2σM=0, C;
    var λ = L, λʹ, iterations = 0;
    do {
        sinλ = Math.sin(λ);
        cosλ = Math.cos(λ);
        sinSqσ = (cosU2*sinλ) * (cosU2*sinλ) + (cosU1*sinU2-sinU1*cosU2*cosλ) * (cosU1*sinU2-sinU1*cosU2*cosλ);
        if (sinSqσ == 0) break;
        sinσ = Math.sqrt(sinSqσ);
        cosσ = sinU1*sinU2 + cosU1*cosU2*cosλ;
        σ = Math.atan2(sinσ, cosσ);
        sinα = cosU1 * cosU2 * sinλ / sinσ;
        cosSqα = 1 - sinα*sinα;
        cos2σM = (cosSqα != 0) ? (cosσ - 2*sinU1*sinU2/cosSqα) : 0;
        C = f/16*cosSqα*(4+f*(4-3*cosSqα));
        λʹ = λ;
        λ = L + (1-C) * f * sinα * (σ + C*sinσ*(cos2σM+C*cosσ*(-1+2*cos2σM*cos2σM)));
        if (Math.abs(λ) > Math.PI) throw new Error('λ > π');
    } while (Math.abs(λ-λʹ) > 1e-12 && ++iterations<1000);
    if (iterations >= 1000) throw new Error('Formula failed to converge');
    var uSq = cosSqα * (a*a - b*b) / (b*b);
    var A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
    var B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));
    var Δσ = B*sinσ*(cos2σM+B/4*(cosσ*(-1+2*cos2σM*cos2σM)-
        B/6*cos2σM*(-3+4*sinσ*sinσ)*(-3+4*cos2σM*cos2σM)));
    var s = b*A*(σ-Δσ);
    var α1 = Math.atan2(cosU2*sinλ,  cosU1*sinU2-sinU1*cosU2*cosλ);
    var α2 = Math.atan2(cosU1*sinλ, -sinU1*cosU2+cosU1*sinU2*cosλ);
    α1 = (α1 + 2*Math.PI) % (2*Math.PI);
    α2 = (α2 + 2*Math.PI) % (2*Math.PI);
    return {
        distance:       s,
        initialBearing: s==0 ? NaN : α1.toDegrees(),
        finalBearing:   s==0 ? NaN : α2.toDegrees(),
        iterations:     iterations,
    };
};
if (Number.prototype.toRadians === undefined) {
    Number.prototype.toRadians = function() { return this * Math.PI / 180; };
}
if (Number.prototype.toDegrees === undefined) {
    Number.prototype.toDegrees = function() { return this * 180 / Math.PI; };
}

function LatLon$3(lat, lon) {
    if (!(this instanceof LatLon$3)) return new LatLon$3(lat, lon);
    this.lat = Number(lat);
    this.lon = Number(lon);
}
LatLon$3.prototype.toVector = function() {
    var φ = this.lat.toRadians();
    var λ = this.lon.toRadians();
    var x = Math.cos(φ) * Math.cos(λ);
    var y = Math.cos(φ) * Math.sin(λ);
    var z = Math.sin(φ);
    return new Vector3d(x, y, z);
};
Vector3d.prototype.toLatLonS = function() {
    var φ = Math.atan2(this.z, Math.sqrt(this.x*this.x + this.y*this.y));
    var λ = Math.atan2(this.y, this.x);
    return new LatLon$3(φ.toDegrees(), λ.toDegrees());
};
LatLon$3.prototype.greatCircle = function(bearing) {
    var φ = this.lat.toRadians();
    var λ = this.lon.toRadians();
    var θ = Number(bearing).toRadians();
    var x =  Math.sin(λ) * Math.cos(θ) - Math.sin(φ) * Math.cos(λ) * Math.sin(θ);
    var y = -Math.cos(λ) * Math.cos(θ) - Math.sin(φ) * Math.sin(λ) * Math.sin(θ);
    var z =  Math.cos(φ) * Math.sin(θ);
    return new Vector3d(x, y, z);
};
Vector3d.prototype.greatCircle = function(bearing) {
    var θ = Number(bearing).toRadians();
    var N = new Vector3d(0, 0, 1);
    var e = N.cross(this);
    var n = this.cross(e);
    var eʹ = e.times(Math.cos(θ)/e.length());
    var nʹ = n.times(Math.sin(θ)/n.length());
    var c = nʹ.minus(eʹ);
    return c;
};
LatLon$3.prototype.distanceTo = function(point, radius) {
    if (!(point instanceof LatLon$3)) throw new TypeError('point is not LatLon object');
    radius = (radius === undefined) ? 6371e3 : Number(radius);
    var p1 = this.toVector();
    var p2 = point.toVector();
    var δ = p1.angleTo(p2);
    var d = δ * radius;
    return d;
};
LatLon$3.prototype.bearingTo = function(point) {
    if (!(point instanceof LatLon$3)) throw new TypeError('point is not LatLon object');
    var p1 = this.toVector();
    var p2 = point.toVector();
    var N = new Vector3d(0, 0, 1);
    var c1 = p1.cross(p2);
    var c2 = p1.cross(N);
    var θ = c1.angleTo(c2, p1);
    return (θ.toDegrees()+360) % 360;
};
LatLon$3.prototype.midpointTo = function(point) {
    if (!(point instanceof LatLon$3)) throw new TypeError('point is not LatLon object');
    var p1 = this.toVector();
    var p2 = point.toVector();
    var mid = p1.plus(p2).unit();
    return mid.toLatLonS();
};
LatLon$3.prototype.intermediatePointTo = function(point, fraction) {
    if (!(point instanceof LatLon$3)) throw new TypeError('point is not LatLon object');
    var n1 = this.toVector();
    var n2 = point.toVector();
    var sinθ = n1.cross(n2).length();
    var cosθ = n1.dot(n2);
    var δ = Math.atan2(sinθ, cosθ);
    var δi = δ * Number(fraction);
    var sinδi = Math.sin(δi);
    var cosδi = Math.cos(δi);
    var d = n1.cross(n2).unit().cross(n1);
    var int = n1.times(cosδi).plus(d.times(sinδi));
    return new Vector3d(int.x, int.y, int.z).toLatLonS();
};
LatLon$3.prototype.intermediatePointOnChordTo = function(point, fraction) {
    if (!(point instanceof LatLon$3)) throw new TypeError('point is not LatLon object');
    var n1 = this.toVector();
    var n2 = point.toVector();
    var int = n1.plus(n2.minus(n1).times(Number(fraction)));
    return new Vector3d(int.x, int.y, int.z).toLatLonS();
};
LatLon$3.prototype.destinationPoint = function(distance, bearing, radius) {
    radius = (radius === undefined) ? 6371e3 : Number(radius);
    var n1 = this.toVector();
    var δ = Number(distance) / radius;
    var θ = Number(bearing).toRadians();
    var N = new Vector3d(0, 0, 1);
    var de = N.cross(n1).unit();
    var dn = n1.cross(de);
    var deSinθ = de.times(Math.sin(θ));
    var dnCosθ = dn.times(Math.cos(θ));
    var d = dnCosθ.plus(deSinθ);
    var x = n1.times(Math.cos(δ));
    var y = d.times(Math.sin(δ));
    var n2 = x.plus(y);
    return n2.toLatLonS();
};
LatLon$3.intersection = function(path1start, path1brngEnd, path2start, path2brngEnd) {
    if (!(path1start instanceof LatLon$3)) throw new TypeError('path1start is not LatLon object');
    if (!(path2start instanceof LatLon$3)) throw new TypeError('path2start is not LatLon object');
    if (!(path1brngEnd instanceof LatLon$3) && isNaN(path1brngEnd)) throw new TypeError('path1brngEnd is not LatLon object or bearing');
    if (!(path2brngEnd instanceof LatLon$3) && isNaN(path2brngEnd)) throw new TypeError('path2brngEnd is not LatLon object or bearing');
    var p1 = path1start.toVector();
    var p2 = path2start.toVector();
    var c1, c2, path1def, path2def;
    if (path1brngEnd instanceof LatLon$3) {
        c1 = p1.cross(path1brngEnd.toVector());
        path1def = 'endpoint';
    } else {
        c1 = path1start.greatCircle(Number(path1brngEnd));
        path1def = 'bearing';
    }
    if (path2brngEnd instanceof LatLon$3) {
        c2 = p2.cross(path2brngEnd.toVector());
        path2def = 'endpoint';
    } else {
        c2 = path2start.greatCircle(Number(path2brngEnd));
        path2def = 'bearing';
    }
    var i1 = c1.cross(c2);
    var i2 = c2.cross(c1);
    var intersection=null, dir1=null, dir2=null;
    switch (path1def+'+'+path2def) {
        case 'bearing+bearing':
            dir1 = Math.sign(c1.cross(p1).dot(i1));
            dir2 = Math.sign(c2.cross(p2).dot(i1));
            switch (dir1+dir2) {
                case  2:
                    intersection = i1;
                    break;
                case -2:
                    intersection = i2;
                    break;
                case  0:
                    intersection = p1.plus(p2).dot(i1) > 0 ? i2 : i1;
                    break;
            }
            break;
        case 'bearing+endpoint':
            dir1 = Math.sign(c1.cross(p1).dot(i1));
            intersection = dir1>0 ? i1 : i2;
            break;
        case 'endpoint+bearing':
            dir2 = Math.sign(c2.cross(p2).dot(i1));
            intersection = dir2>0 ? i1 : i2;
            break;
        case 'endpoint+endpoint':
            var mid = p1.plus(p2).plus(path1brngEnd.toVector()).plus(path2brngEnd.toVector());
            intersection = mid.dot(i1)>0 ? i1 : i2;
            break;
    }
    return intersection.toLatLonS();
};
LatLon$3.prototype.crossTrackDistanceTo = function(pathStart, pathBrngEnd, radius) {
    if (!(pathStart instanceof LatLon$3)) throw new TypeError('pathStart is not LatLon object');
    var R = (radius === undefined) ? 6371e3 : Number(radius);
    var p = this.toVector();
    var gc = pathBrngEnd instanceof LatLon$3
        ? pathStart.toVector().cross(pathBrngEnd.toVector())
        : pathStart.greatCircle(Number(pathBrngEnd));
    var α = gc.angleTo(p) - Math.PI/2;
    var d = α * R;
    return d;
};
LatLon$3.prototype.alongTrackDistanceTo = function(pathStart, pathBrngEnd, radius) {
    if (!(pathStart instanceof LatLon$3)) throw new TypeError('pathStart is not LatLon object');
    var R = (radius === undefined) ? 6371e3 : Number(radius);
    var p = this.toVector();
    var gc = pathBrngEnd instanceof LatLon$3
        ? pathStart.toVector().cross(pathBrngEnd.toVector())
        : pathStart.greatCircle(Number(pathBrngEnd));
    var pat = gc.cross(p).cross(gc);
    var α = pathStart.toVector().angleTo(pat, gc);
    var d = α * R;
    return d;
};
LatLon$3.prototype.nearestPointOnSegment = function(point1, point2) {
    var p = null;
    if (this.isBetween(point1, point2)) {
        var n0 = this.toVector(), n1 = point1.toVector(), n2 = point2.toVector();
        var c1 = n1.cross(n2);
        var c2 = n0.cross(c1);
        var n = c1.cross(c2);
        p = n.toLatLonS();
    } else {
        var d1 = this.distanceTo(point1);
        var d2 = this.distanceTo(point2);
        p = d1<d2 ? point1 : point2;
    }
    return p;
};
LatLon$3.prototype.isBetween = function(point1, point2) {
    var n0 = this.toVector(), n1 = point1.toVector(), n2 = point2.toVector();
    var δ10 = n0.minus(n1), δ12 = n2.minus(n1);
    var δ20 = n0.minus(n2), δ21 = n1.minus(n2);
    var extent1 = δ10.dot(δ12);
    var extent2 = δ20.dot(δ21);
    var isBetween = extent1>=0 && extent2>=0;
    var isSameHemisphere = n0.dot(n1)>=0 && n0.dot(n2)>=0;
    return isBetween && isSameHemisphere;
};
LatLon$3.prototype.enclosedBy = function(polygon) {
    var closed = polygon[0].equals(polygon[polygon.length-1]);
    if (!closed) polygon.push(polygon[0]);
    var nVertices = polygon.length - 1;
    var p = this.toVector();
    var vectorToVertex = [];
    for (var v=0; v<nVertices; v++) vectorToVertex[v] = p.minus(polygon[v].toVector());
    vectorToVertex.push(vectorToVertex[0]);
    var Σθ = 0;
    for (var v=0; v<nVertices; v++) {
        Σθ += vectorToVertex[v].angleTo(vectorToVertex[v+1], p);
    }
    var enclosed = Math.abs(Σθ) > Math.PI;
    if (!closed) polygon.pop();
    return enclosed;
};
LatLon$3.areaOf = function(polygon, radius) {
    var R = (radius == undefined) ? 6371e3 : Number(radius);
    var closed = polygon[0].equals(polygon[polygon.length-1]);
    if (!closed) polygon.push(polygon[0]);
    var n = polygon.length - 1;
    var c = [];
    for (var v=0; v<n; v++) {
        var i = polygon[v].toVector();
        var j = polygon[v+1].toVector();
        c[v] = i.cross(j);
    }
    c.push(c[0]);
    var n1 = polygon[0].toVector();
    var Σα = 0;
    for (var v=0; v<n; v++) Σα += c[v].angleTo(c[v+1], n1);
    var Σθ = n*Math.PI - Math.abs(Σα);
    var E = (Σθ - (n-2)*Math.PI);
    var A = E * R*R;
    if (!closed) polygon.pop();
    return A;
};
LatLon$3.meanOf = function(points) {
    var m = new Vector3d(0, 0, 0);
    for (var p=0; p<points.length; p++) {
        m = m.plus(points[p].toVector());
    }
    return m.unit().toLatLonS();
};
LatLon$3.prototype.equals = function(point) {
    if (!(point instanceof LatLon$3)) throw new TypeError('point is not LatLon object');
    if (this.lat != point.lat) return false;
    if (this.lon != point.lon) return false;
    return true;
};
LatLon$3.prototype.toString = function(format, dp) {
    return Dms.toLat(this.lat, format, dp) + ', ' + Dms.toLon(this.lon, format, dp);
};
if (Number.prototype.toRadians === undefined) {
    Number.prototype.toRadians = function() { return this * Math.PI / 180; };
}
if (Number.prototype.toDegrees === undefined) {
    Number.prototype.toDegrees = function() { return this * 180 / Math.PI; };
}
if (Math.sign === undefined) {
    Math.sign = function(x) {
        x = +x;
        if (x === 0 || isNaN(x)) return x;
        return x > 0 ? 1 : -1;
    };
}

function Utm(zone, hemisphere, easting, northing, datum, convergence, scale) {
    if (!(this instanceof Utm)) {
        return new Utm(zone, hemisphere, easting, northing, datum, convergence, scale);
    }
    if (datum === undefined) datum = LatLon$2.datum.WGS84;
    if (convergence === undefined) convergence = null;
    if (scale === undefined) scale = null;
    if (!(1<=zone && zone<=60)) throw new Error('Invalid UTM zone '+zone);
    if (!hemisphere.match(/[NS]/i)) throw new Error('Invalid UTM hemisphere '+hemisphere);
    this.zone = Number(zone);
    this.hemisphere = hemisphere.toUpperCase();
    this.easting = Number(easting);
    this.northing = Number(northing);
    this.datum = datum;
    this.convergence = convergence===null ? null : Number(convergence);
    this.scale = scale===null ? null : Number(scale);
}
LatLon$2.prototype.toUtm = function() {
    if (isNaN(this.lat) || isNaN(this.lon)) throw new Error('Invalid point');
    if (!(-80<=this.lat && this.lat<=84)) throw new Error('Outside UTM limits');
    var falseEasting = 500e3, falseNorthing = 10000e3;
    var zone = Math.floor((this.lon+180)/6) + 1;
    var λ0 = ((zone-1)*6 - 180 + 3).toRadians();
    var mgrsLatBands = 'CDEFGHJKLMNPQRSTUVWXX';
    var latBand = mgrsLatBands.charAt(Math.floor(this.lat/8+10));
    if (zone==31 && latBand=='V' && this.lon>= 3) { zone++; λ0 += (6).toRadians(); }
    if (zone==32 && latBand=='X' && this.lon<  9) { zone--; λ0 -= (6).toRadians(); }
    if (zone==32 && latBand=='X' && this.lon>= 9) { zone++; λ0 += (6).toRadians(); }
    if (zone==34 && latBand=='X' && this.lon< 21) { zone--; λ0 -= (6).toRadians(); }
    if (zone==34 && latBand=='X' && this.lon>=21) { zone++; λ0 += (6).toRadians(); }
    if (zone==36 && latBand=='X' && this.lon< 33) { zone--; λ0 -= (6).toRadians(); }
    if (zone==36 && latBand=='X' && this.lon>=33) { zone++; λ0 += (6).toRadians(); }
    var φ = this.lat.toRadians();
    var λ = this.lon.toRadians() - λ0;
    var a = this.datum.ellipsoid.a, f = this.datum.ellipsoid.f;
    var k0 = 0.9996;
    var e = Math.sqrt(f*(2-f));
    var n = f / (2 - f);
    var n2 = n*n, n3 = n*n2, n4 = n*n3, n5 = n*n4, n6 = n*n5;
    var cosλ = Math.cos(λ), sinλ = Math.sin(λ), tanλ = Math.tan(λ);
    var τ = Math.tan(φ);
    var σ = Math.sinh(e*Math.atanh(e*τ/Math.sqrt(1+τ*τ)));
    var τʹ = τ*Math.sqrt(1+σ*σ) - σ*Math.sqrt(1+τ*τ);
    var ξʹ = Math.atan2(τʹ, cosλ);
    var ηʹ = Math.asinh(sinλ / Math.sqrt(τʹ*τʹ + cosλ*cosλ));
    var A = a/(1+n) * (1 + 1/4*n2 + 1/64*n4 + 1/256*n6);
    var α = [ null,
        1/2*n - 2/3*n2 + 5/16*n3 +   41/180*n4 -     127/288*n5 +      7891/37800*n6,
              13/48*n2 -  3/5*n3 + 557/1440*n4 +     281/630*n5 - 1983433/1935360*n6,
                       61/240*n3 -  103/140*n4 + 15061/26880*n5 +   167603/181440*n6,
                               49561/161280*n4 -     179/168*n5 + 6601661/7257600*n6,
                                                 34729/80640*n5 - 3418889/1995840*n6,
                                                              212378941/319334400*n6 ];
    var ξ = ξʹ;
    for (var j=1; j<=6; j++) ξ += α[j] * Math.sin(2*j*ξʹ) * Math.cosh(2*j*ηʹ);
    var η = ηʹ;
    for (var j=1; j<=6; j++) η += α[j] * Math.cos(2*j*ξʹ) * Math.sinh(2*j*ηʹ);
    var x = k0 * A * η;
    var y = k0 * A * ξ;
    var pʹ = 1;
    for (var j=1; j<=6; j++) pʹ += 2*j*α[j] * Math.cos(2*j*ξʹ) * Math.cosh(2*j*ηʹ);
    var qʹ = 0;
    for (var j=1; j<=6; j++) qʹ += 2*j*α[j] * Math.sin(2*j*ξʹ) * Math.sinh(2*j*ηʹ);
    var γʹ = Math.atan(τʹ / Math.sqrt(1+τʹ*τʹ)*tanλ);
    var γʺ = Math.atan2(qʹ, pʹ);
    var γ = γʹ + γʺ;
    var sinφ = Math.sin(φ);
    var kʹ = Math.sqrt(1 - e*e*sinφ*sinφ) * Math.sqrt(1 + τ*τ) / Math.sqrt(τʹ*τʹ + cosλ*cosλ);
    var kʺ = A / a * Math.sqrt(pʹ*pʹ + qʹ*qʹ);
    var k = k0 * kʹ * kʺ;
    x = x + falseEasting;
    if (y < 0) y = y + falseNorthing;
    x = Number(x.toFixed(6));
    y = Number(y.toFixed(6));
    var convergence = Number(γ.toDegrees().toFixed(9));
    var scale = Number(k.toFixed(12));
    var h = this.lat>=0 ? 'N' : 'S';
    return new Utm(zone, h, x, y, this.datum, convergence, scale);
};
Utm.prototype.toLatLonE = function() {
    var z = this.zone;
    var h = this.hemisphere;
    var x = this.easting;
    var y = this.northing;
    if (isNaN(z) || isNaN(x) || isNaN(y)) throw new Error('Invalid coordinate');
    var falseEasting = 500e3, falseNorthing = 10000e3;
    var a = this.datum.ellipsoid.a, f = this.datum.ellipsoid.f;
    var k0 = 0.9996;
    x = x - falseEasting;
    y = h=='S' ? y - falseNorthing : y;
    var e = Math.sqrt(f*(2-f));
    var n = f / (2 - f);
    var n2 = n*n, n3 = n*n2, n4 = n*n3, n5 = n*n4, n6 = n*n5;
    var A = a/(1+n) * (1 + 1/4*n2 + 1/64*n4 + 1/256*n6);
    var η = x / (k0*A);
    var ξ = y / (k0*A);
    var β = [ null,
        1/2*n - 2/3*n2 + 37/96*n3 -    1/360*n4 -   81/512*n5 +    96199/604800*n6,
               1/48*n2 +  1/15*n3 - 437/1440*n4 +   46/105*n5 - 1118711/3870720*n6,
                        17/480*n3 -   37/840*n4 - 209/4480*n5 +      5569/90720*n6,
                                 4397/161280*n4 -   11/504*n5 -  830251/7257600*n6,
                                               4583/161280*n5 -  108847/3991680*n6,
                                                             20648693/638668800*n6 ];
    var ξʹ = ξ;
    for (var j=1; j<=6; j++) ξʹ -= β[j] * Math.sin(2*j*ξ) * Math.cosh(2*j*η);
    var ηʹ = η;
    for (var j=1; j<=6; j++) ηʹ -= β[j] * Math.cos(2*j*ξ) * Math.sinh(2*j*η);
    var sinhηʹ = Math.sinh(ηʹ);
    var sinξʹ = Math.sin(ξʹ), cosξʹ = Math.cos(ξʹ);
    var τʹ = sinξʹ / Math.sqrt(sinhηʹ*sinhηʹ + cosξʹ*cosξʹ);
    var τi = τʹ;
    do {
        var σi = Math.sinh(e*Math.atanh(e*τi/Math.sqrt(1+τi*τi)));
        var τiʹ = τi * Math.sqrt(1+σi*σi) - σi * Math.sqrt(1+τi*τi);
        var δτi = (τʹ - τiʹ)/Math.sqrt(1+τiʹ*τiʹ)
            * (1 + (1-e*e)*τi*τi) / ((1-e*e)*Math.sqrt(1+τi*τi));
        τi += δτi;
    } while (Math.abs(δτi) > 1e-12);
    var τ = τi;
    var φ = Math.atan(τ);
    var λ = Math.atan2(sinhηʹ, cosξʹ);
    var p = 1;
    for (var j=1; j<=6; j++) p -= 2*j*β[j] * Math.cos(2*j*ξ) * Math.cosh(2*j*η);
    var q = 0;
    for (var j=1; j<=6; j++) q += 2*j*β[j] * Math.sin(2*j*ξ) * Math.sinh(2*j*η);
    var γʹ = Math.atan(Math.tan(ξʹ) * Math.tanh(ηʹ));
    var γʺ = Math.atan2(q, p);
    var γ = γʹ + γʺ;
    var sinφ = Math.sin(φ);
    var kʹ = Math.sqrt(1 - e*e*sinφ*sinφ) * Math.sqrt(1 + τ*τ) * Math.sqrt(sinhηʹ*sinhηʹ + cosξʹ*cosξʹ);
    var kʺ = A / a / Math.sqrt(p*p + q*q);
    var k = k0 * kʹ * kʺ;
    var λ0 = ((z-1)*6 - 180 + 3).toRadians();
    λ += λ0;
    var lat = Number(φ.toDegrees().toFixed(11));
    var lon = Number(λ.toDegrees().toFixed(11));
    var convergence = Number(γ.toDegrees().toFixed(9));
    var scale = Number(k.toFixed(12));
    var latLong = new LatLon$2(lat, lon, this.datum);
    latLong.convergence = convergence;
    latLong.scale = scale;
    return latLong;
};
Utm.parse = function(utmCoord, datum) {
    if (datum === undefined) datum = LatLon$2.datum.WGS84;
    utmCoord = utmCoord.trim().match(/\S+/g);
    if (utmCoord==null || utmCoord.length!=4) throw new Error('Invalid UTM coordinate ‘'+utmCoord+'’');
    var zone = utmCoord[0], hemisphere = utmCoord[1], easting = utmCoord[2], northing = utmCoord[3];
    return new Utm(zone, hemisphere, easting, northing, datum);
};
Utm.prototype.toString = function(digits) {
    digits = Number(digits||0);
    var z = this.zone<10 ? '0'+this.zone : this.zone;
    var h = this.hemisphere;
    var e = this.easting;
    var n = this.northing;
    if (isNaN(z) || !h.match(/[NS]/) || isNaN(e) || isNaN(n)) return '';
    return z+' '+h+' '+e.toFixed(digits)+' '+n.toFixed(digits);
};
if (Math.sinh === undefined) {
    Math.sinh = function(x) {
        return (Math.exp(x) - Math.exp(-x)) / 2;
    };
}
if (Math.cosh === undefined) {
    Math.cosh = function(x) {
        return (Math.exp(x) + Math.exp(-x)) / 2;
    };
}
if (Math.tanh === undefined) {
    Math.tanh = function(x) {
        return (Math.exp(x) - Math.exp(-x)) / (Math.exp(x) + Math.exp(-x));
    };
}
if (Math.asinh === undefined) {
    Math.asinh = function(x) {
        return Math.log(x + Math.sqrt(1 + x*x));
    };
}
if (Math.atanh === undefined) {
    Math.atanh = function(x) {
        return Math.log((1+x) / (1-x)) / 2;
    };
}

Mgrs.latBands = 'CDEFGHJKLMNPQRSTUVWXX';
Mgrs.e100kLetters = [ 'ABCDEFGH', 'JKLMNPQR', 'STUVWXYZ' ];
Mgrs.n100kLetters = [ 'ABCDEFGHJKLMNPQRSTUV', 'FGHJKLMNPQRSTUVABCDE' ];
function Mgrs(zone, band, e100k, n100k, easting, northing, datum) {
    if (!(this instanceof Mgrs)) return new Mgrs(zone, band, e100k, n100k, easting, northing, datum);
    if (datum === undefined) datum = LatLon$2.datum.WGS84;
    if (!(1<=zone && zone<=60)) throw new Error('Invalid MGRS grid reference (zone ‘'+zone+'’)');
    if (band.length != 1) throw new Error('Invalid MGRS grid reference (band ‘'+band+'’)');
    if (Mgrs.latBands.indexOf(band) == -1) throw new Error('Invalid MGRS grid reference (band ‘'+band+'’)');
    if (e100k.length!=1) throw new Error('Invalid MGRS grid reference (e100k ‘'+e100k+'’)');
    if (n100k.length!=1) throw new Error('Invalid MGRS grid reference (n100k ‘'+n100k+'’)');
    this.zone = Number(zone);
    this.band = band;
    this.e100k = e100k;
    this.n100k = n100k;
    this.easting = Number(easting);
    this.northing = Number(northing);
    this.datum = datum;
}
Utm.prototype.toMgrs = function() {
    if (isNaN(this.zone + this.easting + this.northing)) throw new Error('Invalid UTM coordinate ‘'+this.toString()+'’');
    var zone = this.zone;
    var latlong = this.toLatLonE();
    var band = Mgrs.latBands.charAt(Math.floor(latlong.lat/8+10));
    var col = Math.floor(this.easting / 100e3);
    var e100k = Mgrs.e100kLetters[(zone-1)%3].charAt(col-1);
    var row = Math.floor(this.northing / 100e3) % 20;
    var n100k = Mgrs.n100kLetters[(zone-1)%2].charAt(row);
    var easting = this.easting % 100e3;
    var northing = this.northing % 100e3;
    easting = Number(easting.toFixed(6));
    northing = Number(northing.toFixed(6));
    return new Mgrs(zone, band, e100k, n100k, easting, northing);
};
Mgrs.prototype.toUtm = function() {
    var zone = this.zone;
    var band = this.band;
    var e100k = this.e100k;
    var n100k = this.n100k;
    var easting = this.easting;
    var northing = this.northing;
    var hemisphere = band>='N' ? 'N' : 'S';
    var col = Mgrs.e100kLetters[(zone-1)%3].indexOf(e100k) + 1;
    var e100kNum = col * 100e3;
    var row = Mgrs.n100kLetters[(zone-1)%2].indexOf(n100k);
    var n100kNum = row * 100e3;
    var latBand = (Mgrs.latBands.indexOf(band)-10)*8;
    var nBand = Math.floor(new LatLon$2(latBand, 0).toUtm().northing/100e3)*100e3;
    var n2M = 0;
    while (n2M + n100kNum + northing < nBand) n2M += 2000e3;
    return new Utm(zone, hemisphere, e100kNum+easting, n2M+n100kNum+northing, this.datum);
};
Mgrs.parse = function(mgrsGridRef) {
    mgrsGridRef = mgrsGridRef.trim();
    if (!mgrsGridRef.match(/\s/)) {
        var en = mgrsGridRef.slice(5);
        en = en.slice(0, en.length/2)+' '+en.slice(-en.length/2);
        mgrsGridRef = mgrsGridRef.slice(0, 3)+' '+mgrsGridRef.slice(3, 5)+' '+en;
    }
    mgrsGridRef = mgrsGridRef.match(/\S+/g);
    if (mgrsGridRef==null || mgrsGridRef.length!=4) throw new Error('Invalid MGRS grid reference ‘'+mgrsGridRef+'’');
    var gzd = mgrsGridRef[0];
    var zone = gzd.slice(0, 2);
    var band = gzd.slice(2, 3);
    var en100k = mgrsGridRef[1];
    var e100k = en100k.slice(0, 1);
    var n100k = en100k.slice(1, 2);
    var e = mgrsGridRef[2], n = mgrsGridRef[3];
    e = e.length>=5 ?  e : (e+'00000').slice(0, 5);
    n = n.length>=5 ?  n : (n+'00000').slice(0, 5);
    return new Mgrs(zone, band, e100k, n100k, e, n);
};
Mgrs.prototype.toString = function(digits) {
    digits = (digits === undefined) ? 10 : Number(digits);
    if ([ 2,4,6,8,10 ].indexOf(digits) == -1) throw new Error('Invalid precision ‘'+digits+'’');
    var zone = ('00'+this.zone).slice(-2);
    var band = this.band;
    var e100k = this.e100k;
    var n100k = this.n100k;
    var eRounded = Math.floor(this.easting/Math.pow(10, 5-digits/2));
    var nRounded = Math.floor(this.northing/Math.pow(10, 5-digits/2));
    var easting = ('00000'+eRounded).slice(-digits/2);
    var northing = ('00000'+nRounded).slice(-digits/2);
    return zone+band + ' ' + e100k+n100k + ' '  + easting + ' ' + northing;
};

function OsGridRef(easting, northing) {
    if (!(this instanceof OsGridRef)) return new OsGridRef(easting, northing);
    this.easting = Number(easting);
    this.northing = Number(northing);
}
OsGridRef.latLonToOsGrid = function(point) {
    if (!(point instanceof LatLon$2)) throw new TypeError('point is not LatLon object');
    if (point.datum != LatLon$2.datum.OSGB36) point = point.convertDatum(LatLon$2.datum.OSGB36);
    var φ = point.lat.toRadians();
    var λ = point.lon.toRadians();
    var a = 6377563.396, b = 6356256.909;
    var F0 = 0.9996012717;
    var φ0 = (49).toRadians(), λ0 = (-2).toRadians();
    var N0 = -100000, E0 = 400000;
    var e2 = 1 - (b*b)/(a*a);
    var n = (a-b)/(a+b), n2 = n*n, n3 = n*n*n;
    var cosφ = Math.cos(φ), sinφ = Math.sin(φ);
    var ν = a*F0/Math.sqrt(1-e2*sinφ*sinφ);
    var ρ = a*F0*(1-e2)/Math.pow(1-e2*sinφ*sinφ, 1.5);
    var η2 = ν/ρ-1;
    var Ma = (1 + n + (5/4)*n2 + (5/4)*n3) * (φ-φ0);
    var Mb = (3*n + 3*n*n + (21/8)*n3) * Math.sin(φ-φ0) * Math.cos(φ+φ0);
    var Mc = ((15/8)*n2 + (15/8)*n3) * Math.sin(2*(φ-φ0)) * Math.cos(2*(φ+φ0));
    var Md = (35/24)*n3 * Math.sin(3*(φ-φ0)) * Math.cos(3*(φ+φ0));
    var M = b * F0 * (Ma - Mb + Mc - Md);
    var cos3φ = cosφ*cosφ*cosφ;
    var cos5φ = cos3φ*cosφ*cosφ;
    var tan2φ = Math.tan(φ)*Math.tan(φ);
    var tan4φ = tan2φ*tan2φ;
    var I = M + N0;
    var II = (ν/2)*sinφ*cosφ;
    var III = (ν/24)*sinφ*cos3φ*(5-tan2φ+9*η2);
    var IIIA = (ν/720)*sinφ*cos5φ*(61-58*tan2φ+tan4φ);
    var IV = ν*cosφ;
    var V = (ν/6)*cos3φ*(ν/ρ-tan2φ);
    var VI = (ν/120) * cos5φ * (5 - 18*tan2φ + tan4φ + 14*η2 - 58*tan2φ*η2);
    var Δλ = λ-λ0;
    var Δλ2 = Δλ*Δλ, Δλ3 = Δλ2*Δλ, Δλ4 = Δλ3*Δλ, Δλ5 = Δλ4*Δλ, Δλ6 = Δλ5*Δλ;
    var N = I + II*Δλ2 + III*Δλ4 + IIIA*Δλ6;
    var E = E0 + IV*Δλ + V*Δλ3 + VI*Δλ5;
    N = Number(N.toFixed(3));
    E = Number(E.toFixed(3));
    return new OsGridRef(E, N);
};
OsGridRef.osGridToLatLon = function(gridref, datum) {
    if (!(gridref instanceof OsGridRef)) throw new TypeError('gridref is not OsGridRef object');
    if (datum === undefined) datum = LatLon$2.datum.WGS84;
    var E = gridref.easting;
    var N = gridref.northing;
    var a = 6377563.396, b = 6356256.909;
    var F0 = 0.9996012717;
    var φ0 = (49).toRadians(), λ0 = (-2).toRadians();
    var N0 = -100000, E0 = 400000;
    var e2 = 1 - (b*b)/(a*a);
    var n = (a-b)/(a+b), n2 = n*n, n3 = n*n*n;
    var φ=φ0, M=0;
    do {
        φ = (N-N0-M)/(a*F0) + φ;
        var Ma = (1 + n + (5/4)*n2 + (5/4)*n3) * (φ-φ0);
        var Mb = (3*n + 3*n*n + (21/8)*n3) * Math.sin(φ-φ0) * Math.cos(φ+φ0);
        var Mc = ((15/8)*n2 + (15/8)*n3) * Math.sin(2*(φ-φ0)) * Math.cos(2*(φ+φ0));
        var Md = (35/24)*n3 * Math.sin(3*(φ-φ0)) * Math.cos(3*(φ+φ0));
        M = b * F0 * (Ma - Mb + Mc - Md);
    } while (N-N0-M >= 0.00001);
    var cosφ = Math.cos(φ), sinφ = Math.sin(φ);
    var ν = a*F0/Math.sqrt(1-e2*sinφ*sinφ);
    var ρ = a*F0*(1-e2)/Math.pow(1-e2*sinφ*sinφ, 1.5);
    var η2 = ν/ρ-1;
    var tanφ = Math.tan(φ);
    var tan2φ = tanφ*tanφ, tan4φ = tan2φ*tan2φ, tan6φ = tan4φ*tan2φ;
    var secφ = 1/cosφ;
    var ν3 = ν*ν*ν, ν5 = ν3*ν*ν, ν7 = ν5*ν*ν;
    var VII = tanφ/(2*ρ*ν);
    var VIII = tanφ/(24*ρ*ν3)*(5+3*tan2φ+η2-9*tan2φ*η2);
    var IX = tanφ/(720*ρ*ν5)*(61+90*tan2φ+45*tan4φ);
    var X = secφ/ν;
    var XI = secφ/(6*ν3)*(ν/ρ+2*tan2φ);
    var XII = secφ/(120*ν5)*(5+28*tan2φ+24*tan4φ);
    var XIIA = secφ/(5040*ν7)*(61+662*tan2φ+1320*tan4φ+720*tan6φ);
    var dE = (E-E0), dE2 = dE*dE, dE3 = dE2*dE, dE4 = dE2*dE2, dE5 = dE3*dE2, dE6 = dE4*dE2, dE7 = dE5*dE2;
    φ = φ - VII*dE2 + VIII*dE4 - IX*dE6;
    var λ = λ0 + X*dE - XI*dE3 + XII*dE5 - XIIA*dE7;
    var point =  new LatLon$2(φ.toDegrees(), λ.toDegrees(), LatLon$2.datum.OSGB36);
    if (datum != LatLon$2.datum.OSGB36) point = point.convertDatum(datum);
    return point;
};
OsGridRef.parse = function(gridref) {
    gridref = String(gridref).trim();
    var match = gridref.match(/^(\d+),\s*(\d+)$/);
    if (match) return new OsGridRef(match[1], match[2]);
    match = gridref.match(/^[A-Z]{2}\s*[0-9]+\s*[0-9]+$/i);
    if (!match) throw new Error('Invalid grid reference');
    var l1 = gridref.toUpperCase().charCodeAt(0) - 'A'.charCodeAt(0);
    var l2 = gridref.toUpperCase().charCodeAt(1) - 'A'.charCodeAt(0);
    if (l1 > 7) l1--;
    if (l2 > 7) l2--;
    var e100km = ((l1-2)%5)*5 + (l2%5);
    var n100km = (19-Math.floor(l1/5)*5) - Math.floor(l2/5);
    var en = gridref.slice(2).trim().split(/\s+/);
    if (en.length == 1) en = [ en[0].slice(0, en[0].length/2), en[0].slice(en[0].length/2) ];
    if (e100km<0 || e100km>6 || n100km<0 || n100km>12) throw new Error('Invalid grid reference');
    if (en.length != 2) throw new Error('Invalid grid reference');
    if (en[0].length != en[1].length) throw new Error('Invalid grid reference');
    en[0] = (en[0]+'00000').slice(0, 5);
    en[1] = (en[1]+'00000').slice(0, 5);
    var e = e100km + en[0];
    var n = n100km + en[1];
    return new OsGridRef(e, n);
};
OsGridRef.prototype.toString = function(digits) {
    digits = (digits === undefined) ? 10 : Number(digits);
    if (isNaN(digits) || digits%2!=0 || digits>16) throw new RangeError('Invalid precision ‘'+digits+'’');
    var e = this.easting;
    var n = this.northing;
    if (isNaN(e) || isNaN(n)) throw new Error('Invalid grid reference');
    if (digits == 0) {
        var eInt = Math.floor(e), eDec = e - eInt;
        var nInt = Math.floor(n), nDec = n - nInt;
        var ePad = ('000000'+eInt).slice(-6) + (eDec>0 ? eDec.toFixed(3).slice(1) : '');
        var nPad = (nInt<1e6 ? ('000000'+nInt).slice(-6) : nInt) + (nDec>0 ? nDec.toFixed(3).slice(1) : '');
        return ePad + ',' + nPad;
    }
    var e100k = Math.floor(e/100000), n100k = Math.floor(n/100000);
    if (e100k<0 || e100k>6 || n100k<0 || n100k>12) return '';
    var l1 = (19-n100k) - (19-n100k)%5 + Math.floor((e100k+10)/5);
    var l2 = (19-n100k)*5%25 + e100k%5;
    if (l1 > 7) l1++;
    if (l2 > 7) l2++;
    var letterPair = String.fromCharCode(l1+'A'.charCodeAt(0), l2+'A'.charCodeAt(0));
    e = Math.floor((e%100000)/Math.pow(10, 5-digits/2));
    n = Math.floor((n%100000)/Math.pow(10, 5-digits/2));
    e = ('00000000'+e).slice(-digits/2);
    n = ('00000000'+n).slice(-digits/2);
    return letterPair + ' ' + e + ' ' + n;
};

for (const prop in LatLon$2) {
  LatLon$2[prop] = LatLon$2[prop];
}

export { LatLon as LatLonSpherical, LatLon$2 as LatLonEllipsoidal, LatLon$3 as LatLonVectors, Vector3d, Utm, Mgrs, OsGridRef, Dms };
