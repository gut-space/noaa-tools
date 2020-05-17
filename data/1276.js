
    // Code generated with noaa-tools.
    // This can be executed using Cesium. The easiest way to test it is to go to
    // https://sandcastle.cesium.com and paste this code there.

    var viewer = new Cesium.Viewer('cesiumContainer');
    var pinBuilder = new Cesium.PinBuilder();
     var czml = [
{
    "id": "document",
    "version": "1.0",
    "name": "document_packet",
    "clock": {
        "interval": "2020-04-12T08:56:03Z/2020-04-12T09:12:06Z",
        "currentTime": "2020-04-12T08:56:03Z",
        "multiplier": 60,
        "range": "LOOP_STOP",
        "step": "SYSTEM_CLOCK_MULTIPLIER"
    }
},
{
    "id": "custom_properties",
    "properties": {
        "custom_attractor": true,
        "ellipsoid": [
            {
                "array": [
                    6378136.6,
                    6378136.6,
                    6356751.9
                ]
            }
        ],
        "map_url": [
            "https://upload.wikimedia.org/wikipedia/commons/c/c4/Earthmap1000x500compac.jpg"
        ],
        "scene3D": true
    }
},
{
    "id": 0,
    "availability": "2020-04-12T08:56:03Z/2020-04-12T09:12:06Z",
    "position": {
        "epoch": "2020-04-12T08:56:03Z",
        "interpolationAlgorithm": "LAGRANGE",
        "interpolationDegree": 5,
        "referenceFrame": "INERTIAL",
        "cartesian": [
            0.0,
            546809.6949856463,
            997167.3135742369,
            7128330.303366815,
            96.3403478,
            1188103.644765849,
            681079.352082901,
            7087761.564471393,
            192.6806956,
            1817723.0522245967,
            358298.7724538702,
            6977545.958026942,
            289.0210434,
            2429484.5913164364,
            31997.997493082,
            6798779.905981012,
            385.3613912,
            3017383.697285734,
            -294617.2016414291,
            6553235.350300109,
            481.701739,
            3575654.0392880384,
            -618339.5775687741,
            6243341.522087208,
            578.0420868,
            4098824.475112825,
            -935992.2298880289,
            5872160.12869854,
            674.3824346,
            4581772.910096309,
            -1244460.0508322448,
            5443354.242960914,
            770.7227824,
            5019776.520857897,
            -1540720.471675621,
            4961151.244429712,
            867.0631302,
            5408557.848776842,
            -1821873.1988128214,
            4430300.223215314,
            963.403478,
            5744326.317360932,
            -2085168.642480352,
            3856024.311564967,
            1059.7438258,
            6023814.781011917,
            -2328034.7582685337,
            3243968.4565759734
        ]
    },
    "billboard": {
        "image": "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAAAXNSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAADJSURBVDhPnZHRDcMgEEMZjVEYpaNklIzSEfLfD4qNnXAJSFWfhO7w2Zc0Tf9QG2rXrEzSUeZLOGm47WoH95x3Hl3jEgilvDgsOQUTqsNl68ezEwn1vae6lceSEEYvvWNT/Rxc4CXQNGadho1NXoJ+9iaqc2xi2xbt23PJCDIB6TQjOC6Bho/sDy3fBQT8PrVhibU7yBFcEPaRxOoeTwbwByCOYf9VGp1BYI1BA+EeHhmfzKbBoJEQwn1yzUZtyspIQUha85MpkNIXB7GizqDEECsAAAAASUVORK5CYII=",
        "show": true
    },
    "label": {
        "text": "NOAA 18",
        "font": "11pt Lucida Console",
        "style": "FILL",
        "fillColor": {
            "rgba": [
                255,
                255,
                0,
                255
            ]
        },
        "outlineColor": {
            "rgba": [
                255,
                255,
                0,
                255
            ]
        },
        "outlineWidth": 1.0
    },
    "path": {
        "show": true,
        "width": 3,
        "resolution": 120,
        "material": {
            "solidColor": {
                "color": {
                    "rgba": [
                        125,
                        80,
                        120,
                        255
                    ]
                }
            }
        }
    }
},
];
var dataSourcePromise = viewer.dataSources.add(Cesium.CzmlDataSource.load(czml));


    var questionPin = viewer.entities.add({
        name : 'AOS(PYMAP3D) lat=79.999313,lon=66.591512',
        position : Cesium.Cartesian3.fromDegrees(66.591512, 79.999313, 0.000000),
        billboard : {
            image : pinBuilder.fromText('A', Cesium.Color.BLUE, 48).toDataURL(),
            verticalOrigin : Cesium.VerticalOrigin.BOTTOM
        }
    });
    
    var questionPin = viewer.entities.add({
        name : 'LOS(PYMAP3D) lat=28.471957,lon=4.010928',
        position : Cesium.Cartesian3.fromDegrees(4.010928, 28.471957, 0.000000),
        billboard : {
            image : pinBuilder.fromText('L', Cesium.Color.YELLOW, 48).toDataURL(),
            verticalOrigin : Cesium.VerticalOrigin.BOTTOM
        }
    });
    
    var questionPin = viewer.entities.add({
        name : 'Upper Left lat=84.234023,lon=43.476290',
        position : Cesium.Cartesian3.fromDegrees(43.476290, 84.234023, 0.000000),
        billboard : {
            image : pinBuilder.fromText('U', Cesium.Color.RED, 48).toDataURL(),
            verticalOrigin : Cesium.VerticalOrigin.BOTTOM
        }
    });
    
    var questionPin = viewer.entities.add({
        name : 'Upper Right lat=75.135412,lon=75.435497',
        position : Cesium.Cartesian3.fromDegrees(75.435497, 75.135412, 0.000000),
        billboard : {
            image : pinBuilder.fromText('U', Cesium.Color.GREEN, 48).toDataURL(),
            verticalOrigin : Cesium.VerticalOrigin.BOTTOM
        }
    });
    
    var questionPin = viewer.entities.add({
        name : 'Lower Left lat=29.583794,lon=-1.810557',
        position : Cesium.Cartesian3.fromDegrees(-1.810557, 29.583794, 0.000000),
        billboard : {
            image : pinBuilder.fromText('L', Cesium.Color.RED, 48).toDataURL(),
            verticalOrigin : Cesium.VerticalOrigin.BOTTOM
        }
    });
    
    var questionPin = viewer.entities.add({
        name : 'Lower Right lat=27.117883,lon=9.698307',
        position : Cesium.Cartesian3.fromDegrees(9.698307, 27.117883, 0.000000),
        billboard : {
            image : pinBuilder.fromText('L', Cesium.Color.GREEN, 48).toDataURL(),
            verticalOrigin : Cesium.VerticalOrigin.BOTTOM
        }
    });
    