
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
        "interval": "2020-04-12T09:01:03.063/2020-04-12T09:17:06.467",
        "currentTime": "2020-04-12T09:01:03",
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
                    6378137.0,
                    6378137.0,
                    6356752.314245179
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
    "availability": "2020-04-12T09:01:03/2020-04-12T09:17:06",
    "position": {
        "epoch": "2020-04-12T09:01:03.063",
        "interpolationAlgorithm": "LAGRANGE",
        "interpolationDegree": 5,
        "referenceFrame": "INERTIAL",
        "cartesian": [
            0.0,
            2497796.7943112236,
            -5267.4637026207,
            6774129.2300034035,
            96.3403478,
            3082605.484606515,
            -331714.630971391,
            6521116.680282644,
            192.6806956,
            3637145.9767669328,
            -654904.8089422454,
            6204072.367708097,
            289.0210434,
            4155984.2504608897,
            -971666.5618236214,
            5826128.216822845,
            385.3613912,
            4634041.249626315,
            -1278893.791563381,
            5391013.587532978,
            481.701739,
            5066642.6596350605,
            -1573576.3845902157,
            4903017.5635137595,
            578.0420868,
            5449564.500818173,
            -1852829.8136067307,
            4366945.753554468,
            674.3824346,
            5779074.098392071,
            -2113923.3992100493,
            3788072.0768559477,
            770.7227824,
            6051966.042593345,
            -2354306.953629566,
            3172086.050718391,
            867.0631302,
            6265592.809754187,
            -2571635.5491272644,
            2525036.13957301,
            963.403478,
            6417889.774265634,
            -2763792.176230177,
            1853269.757675314,
            1059.7438258,
            6507394.402034201,
            -2928908.081567771,
            1163370.5438812037
        ]
    },
    "billboard": {
        "image": "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAAAXNSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAADJSURBVDhPnZHRDcMgEEMZjVEYpaNklIzSEfLfD4qNnXAJSFWfhO7w2Zc0Tf9QG2rXrEzSUeZLOGm47WoH95x3Hl3jEgilvDgsOQUTqsNl68ezEwn1vae6lceSEEYvvWNT/Rxc4CXQNGadho1NXoJ+9iaqc2xi2xbt23PJCDIB6TQjOC6Bho/sDy3fBQT8PrVhibU7yBFcEPaRxOoeTwbwByCOYf9VGp1BYI1BA+EeHhmfzKbBoJEQwn1yzUZtyspIQUha85MpkNIXB7GizqDEECsAAAAASUVORK5CYII=",
        "show": true
    },
    "label": {
        "text": "satname",
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

 // points 0 out of 3
    var questionPin = viewer.entities.add({
        name : 'AOS:2020-04-12 09:01:03.063476+00:00, method spherical',
        position : Cesium.Cartesian3.fromDegrees(79.859318, 40.585883, 0.000000),
        billboard : {
            image : pinBuilder.fromText('A', Cesium.Color.RED, 48).toDataURL(),
            verticalOrigin : Cesium.VerticalOrigin.BOTTOM
        }
    });
    
    var questionPin = viewer.entities.add({
        name : 'LOS:2020-04-12 09:17:06.466954+00:00, method spherical',
        position : Cesium.Cartesian3.fromDegrees(28.223308, -18.162042, 0.000000),
        billboard : {
            image : pinBuilder.fromText('L', Cesium.Color.RED, 48).toDataURL(),
            verticalOrigin : Cesium.VerticalOrigin.BOTTOM
        }
    });
     // points 1 out of 3
    var questionPin = viewer.entities.add({
        name : 'AOS:2020-04-12 09:01:03.063476+00:00, method oblate',
        position : Cesium.Cartesian3.fromDegrees(79.917969, 40.585883, 0.000000),
        billboard : {
            image : pinBuilder.fromText('A', Cesium.Color.RED, 48).toDataURL(),
            verticalOrigin : Cesium.VerticalOrigin.BOTTOM
        }
    });
    
    var questionPin = viewer.entities.add({
        name : 'LOS:2020-04-12 09:17:06.466954+00:00, method oblate',
        position : Cesium.Cartesian3.fromDegrees(28.364909, -18.162042, 0.000000),
        billboard : {
            image : pinBuilder.fromText('L', Cesium.Color.RED, 48).toDataURL(),
            verticalOrigin : Cesium.VerticalOrigin.BOTTOM
        }
    });
     // points 2 out of 3
    var questionPin = viewer.entities.add({
        name : 'AOS:2020-04-12 09:01:03.063476+00:00, method pymap3d',
        position : Cesium.Cartesian3.fromDegrees(79.999313, 66.591512, 0.000000),
        billboard : {
            image : pinBuilder.fromText('A', Cesium.Color.RED, 48).toDataURL(),
            verticalOrigin : Cesium.VerticalOrigin.BOTTOM
        }
    });
    
    var questionPin = viewer.entities.add({
        name : 'LOS:2020-04-12 09:17:06.466954+00:00, method pymap3d',
        position : Cesium.Cartesian3.fromDegrees(28.471957, 4.010928, 0.000000),
        billboard : {
            image : pinBuilder.fromText('L', Cesium.Color.RED, 48).toDataURL(),
            verticalOrigin : Cesium.VerticalOrigin.BOTTOM
        }
    });
    