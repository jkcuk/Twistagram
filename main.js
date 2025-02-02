
// This code is based on three.js, which comes with the following license:
//
// The MIT License
//
// Copyright © 2010-2024 three.js authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
import * as THREE from 'three';

import { GUI } from 'three/addons/libs/lil-gui.module.min.js';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

import { DragControls } from 'three/addons/controls/DragControls.js';
import { VRButton } from 'three/addons/webxr/VRButton.js';
import { HTMLMesh } from 'three/addons/interactive/HTMLMesh.js';
import { InteractiveGroup } from 'three/addons/interactive/InteractiveGroup.js';
import { XRControllerModelFactory } from 'three/addons/webxr/XRControllerModelFactory.js';

let name = 'Twistagram';

// variables that are not uniforms
let aspectRatioVideoFeedU = 4.0/3.0;
let aspectRatioVideoFeedE = 4.0/3.0;

// Nokia HR20, according to https://www.camerafv5.com/devices/manufacturers/hmd_global/nokia_xr20_ttg_0/
let fovVideoFeedU = 67.3;	// (user-facing) camera
let fovVideoFeedE = 68.3;	// (environment-facing) camera
// the FOV of the screen depends on the user's distance from the screen, of course
let fovScreen = 140;	// approximation to the horizontal visual field -- see https://en.wikipedia.org/wiki/Visual_field

let ipd = 0.064;	// interpupillary distance, see https://en.wikipedia.org/wiki/Pupillary_distance#Databases
let componentDistance = 0.2;
let alpha = 0.0;
let centreOfObjectPlane = new THREE.Vector3(0, 0, -10);	// in component's coordinate system
let designViewPosition = new THREE.Vector3(0, 0, 0.04);	// in component's coordinate system
let inFrontOfCamera = false;

// camera with wide aperture
let focusDistance = 1e8;
let noOfRays = 1;

let raytracingSphereRadius = 100.0;	// what this code does should be independent of this value, but isn't (?)

// "internal" variables
let raytracingSphere;
let raytracingSphereShaderMaterial;

let scene;
let renderer;
let videoFeedU, videoFeedE;	// feeds from user/environment-facing cameras
let camera;
let orbitControls;
let dragControls;
let background = 1;
let backgroundTexture;

// the menu
let gui;
let GUIMesh;
let visibleControl, perfectRotatorVisibleControl, vrControlsVisibleControl, backgroundControl, componentYControl, inFrontOfCameraControl;

// lift the component up to eye level (in case of VR only)
let componentY = 0.0;

// the status text area
let status = document.createElement('div');
let statusTime;	// the time the last status was posted

// the info text area
let info = document.createElement('div');

// the stored photo
let storedPhoto;
let storedPhotoDescription;
let storedPhotoInfoString;
let showingStoredPhoto = false;

// my Canon EOS450D, recorded with the iPad, sent to the Mac, and then exported as Audio only
const click = new Audio('./click.m4a');

// a piece of paper being crumpled, recorded with the iPad, sent to the Mac, and then exported as Audio only
const bin = new Audio('./bin.m4a');

init();
animate();

function init() {
	// create the info element first so that any problems can be communicated
	createStatus();

	scene = new THREE.Scene();
	// scene.background = new THREE.Color( 'skyblue' );
	let windowAspectRatio = window.innerWidth / window.innerHeight;
	camera = new THREE.PerspectiveCamera( fovScreen, windowAspectRatio, 0.01, 2*raytracingSphereRadius );
	camera.position.z = componentDistance;
	screenChanged();
	
	renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
	renderer.setPixelRatio(window.devicePixelRatio);
	renderer.setSize( window.innerWidth, window.innerHeight );
	// document.body.appendChild( renderer.domElement );
	renderer.xr.enabled = true;
	document.body.appendChild( VRButton.createButton( renderer ) );	// for VR content
	document.body.appendChild( renderer.domElement );
	// document.getElementById('livePhoto').appendChild( renderer.domElement );

	loadBackgroundImage();

	createVideoFeeds();

	addRaytracingSphere();

	// user interface

	addEventListenersEtc();

	addOrbitControls();

	// the controls menu
	createGUI();

	addDragControls();

	// check if VR is supported (see https://developer.mozilla.org/en-US/docs/Web/API/XRSystem/isSessionSupported)...
	// if (navigator.xr) {
	if ( 'xr' in navigator ) {
		// renderer.xr.enabled = false;
		// navigator.xr.isSessionSupported("immersive-vr").then((isSupported) => {
		navigator.xr.isSessionSupported( 'immersive-vr' ).then( function ( supported ) {
			if (supported) {
				// ... and enable the relevant features
				renderer.xr.enabled = true;
				// use renderer.xr.isPresenting to find out if we are in XR mode -- see https://threejs.org/docs/#api/en/renderers/webxr/WebXRManager 
				// (and https://threejs.org/docs/#api/en/renderers/WebGLRenderer.xr, which states that renderer.xr points to the WebXRManager)
				document.body.appendChild( VRButton.createButton( renderer ) );	// for VR content
				addXRInteractivity();
			}
		});
	}

	createInfo();
	refreshInfo();
}

function animate() {
	renderer.setAnimationLoop( render );
}

function render() {
	// requestAnimationFrame( animate );

	// stats.begin();

	if(!showingStoredPhoto) {
		// update uniforms
		updateUniforms();

		// ensure the raytracing  sphere is centred on  the camera
		raytracingSphere.position.copy(camera.position);

		renderer.render( scene,  camera );
	}

	// stats.end();
}

function updateUniforms() {
	// // the component is centred at (uniform) centreOfComponent, with basis vectors aHat, bHat, cHat;
	// // create a 4-matrix that converts from (a, b, c) to (x, y, z)
	// let q = new THREE.Quaternion();
	// q.setFromAxisAngle( axis, angleRad );

	// let m = new THREE.Matrix4();
	// m.compose(
	// 	centreOfComponent,	// position
	// 	q,	// rotation
	// 	new THREE.Vector3(1, 1, 1)	// scale
	// )
	// m.makeBasis(aHat, bHat, cHat);

	// the tangents for the environment-facing camera video feed
	let tanHalfFovHE, tanHalfFovVE;
	if(aspectRatioVideoFeedE > 1.0) {
		// horizontal orientation
		tanHalfFovHE = Math.tan(0.5*fovVideoFeedE*Math.PI/180.0);
		tanHalfFovVE = Math.tan(0.5*fovVideoFeedE*Math.PI/180.0)/aspectRatioVideoFeedE;
	} else {
		// vertical orientation
		tanHalfFovHE = Math.tan(0.5*fovVideoFeedE*Math.PI/180.0)*aspectRatioVideoFeedE;
		tanHalfFovVE = Math.tan(0.5*fovVideoFeedE*Math.PI/180.0);
	}
	raytracingSphereShaderMaterial.uniforms.halfWidthE.value = raytracingSphereShaderMaterial.uniforms.videoDistance.value*tanHalfFovHE;
	raytracingSphereShaderMaterial.uniforms.halfHeightE.value = raytracingSphereShaderMaterial.uniforms.videoDistance.value*tanHalfFovVE;

	// the tangents for the user-facing camera video feed
	let tanHalfFovHU, tanHalfFovVU;
	if(aspectRatioVideoFeedU > 1.0) {
		// horizontal orientation
		tanHalfFovHU = Math.tan(0.5*fovVideoFeedU*Math.PI/180.0);
		tanHalfFovVU = Math.tan(0.5*fovVideoFeedU*Math.PI/180.0)/aspectRatioVideoFeedU;
	} else {
		// vertical orientation
		tanHalfFovHU = Math.tan(0.5*fovVideoFeedU*Math.PI/180.0)*aspectRatioVideoFeedU;
		tanHalfFovVU = Math.tan(0.5*fovVideoFeedU*Math.PI/180.0);
	}
	raytracingSphereShaderMaterial.uniforms.halfWidthU.value = raytracingSphereShaderMaterial.uniforms.videoDistance.value*tanHalfFovHU;
	raytracingSphereShaderMaterial.uniforms.halfHeightU.value = raytracingSphereShaderMaterial.uniforms.videoDistance.value*tanHalfFovVU;

	raytracingSphereShaderMaterial.uniforms.showVideoFeed.value = (background == 0);
	raytracingSphereShaderMaterial.uniforms.backgroundTexture.value = backgroundTexture;

	// the camera's local coordinate system
	let cameraXHat = new THREE.Vector3();	// right direction
	let cameraYHat = new THREE.Vector3();	// up direction
	let cameraZHat = new THREE.Vector3();	// backwards direction, i.e. - viewDirection
	camera.matrixWorld.extractBasis(cameraXHat, cameraYHat, cameraZHat);

/* 	// alternative way:
	let viewDirection = new THREE.Vector3();
	camera.getWorldDirection(viewDirection);
	// the standard way of establishing the "up" direction is to calculate xHat = viewDirection x (0, 1, 0)
	// and then calculating yHat = xHat x viewDirection;
	// we need to deal with the special case when viewDirection \propto (0, 1, 0)
	if((viewDirection.x == 0.0) && (viewDirection.z == 0.0)) {
		// viewDirection is along y direction; make x the x direction, and then the y direction is +z or -z
		cameraXHat = new THREE.Vector3(1, 0, 0);
	} else {
		// viewDirection is not along y direction
		cameraXHat.crossVectors(viewDirection, new THREE.Vector3(0, 1, 0)).normalize();
	}
	cameraYHat.crossVectors(cameraXHat, viewDirection).normalize();
 */

	// calculate the component's model matrix and model-matrix inverse

	// calculate the component's position
	let componentPosition = new THREE.Vector3();	
	// let comparisonPosition = new THREE.Vector3();	
	componentPosition.copy(camera.position);	// the camera position
	componentPosition.addScaledVector(cameraZHat, -componentDistance);	// -componentDistance * zHat
	componentPosition.addScaledVector(cameraXHat, 0.5*ipd);	// a bit to the right
	// comparisonPosition.copy(componentPosition);
	// comparisonPosition.addScaledVector(cameraXHat, -0.5*ipd);	// a bit to the left

	// calculate the component's model matrix and its inverse
	let componentModelMatrix = new THREE.Matrix4();
	let componentModelMatrixInverse = new THREE.Matrix4();
	// let comparisonModelMatrix = new THREE.Matrix4();
	// let comparisonModelMatrixInverse = new THREE.Matrix4();
	if(inFrontOfCamera) {
		componentModelMatrix.copy(camera.matrixWorld);	// the basis vectors are the same as those in the camera's model matrix...
		componentModelMatrix.setPosition(componentPosition);	// ... but the position isn't the same

		// comparisonModelMatrix.copy(camera.matrixWorld);	// the basis vectors are the same as those in the camera's model matrix...
		// comparisonModelMatrix.setPosition(comparisonPosition);	// ... but the position isn't the same
	
		// comparisonModelMatrixInverse.copy(comparisonModelMatrix);	// start from the model matrix...
		// comparisonModelMatrixInverse.invert();	// ... and invert it
	} else {
		componentModelMatrix.makeTranslation(0, componentY, 0);

		// comparisonModelMatrix.identity();
		// comparisonModelMatrixInverse.identity();
	}
		
	componentModelMatrixInverse.copy(componentModelMatrix);	// start from the model matrix...
	componentModelMatrixInverse.invert();	// ... and invert it

	raytracingSphereShaderMaterial.uniforms.componentModelMatrix.value = componentModelMatrix;	// set the uniform
	raytracingSphereShaderMaterial.uniforms.componentModelMatrixInverse.value = componentModelMatrixInverse;	// set the uniform
	// raytracingSphereShaderMaterial.uniforms.comparisonModelMatrix.value = comparisonModelMatrix;	// set the uniform
	// raytracingSphereShaderMaterial.uniforms.comparisonModelMatrixInverse.value = comparisonModelMatrixInverse;	// set the uniform

	// if(counter < 10) console.log(`viewDirection = (${viewDirection.x.toPrecision(2)}, ${viewDirection.y.toPrecision(2)}, ${viewDirection.z.toPrecision(2)})`);
	// raytracingSphereShaderMaterial.uniforms.viewDirection.value = cameraZHat.multiplyScalar(-1);

	// if the component is "glued" to the eye, place it in front of the camera
	// raytracingSphereShaderMaterial.uniforms.centreOfComponent.value.y = componentY;	// camera.position.y;
	// raytracingSphereShaderMaterial.uniforms.designViewPosition.value.y = componentY;
	// raytracingSphereShaderMaterial.uniforms.centreOfPerfectRotator.value.y = componentY;
	// raytracingSphereShaderMaterial.uniforms.centreOfPerfectRotator.value.x = componentPosition.x - 0.064;


	raytracingSphereShaderMaterial.uniforms.noOfRays.value = noOfRays;
	raytracingSphereShaderMaterial.uniforms.cameraXHat.value.copy(cameraXHat);
	raytracingSphereShaderMaterial.uniforms.cameraYHat.value.copy(cameraYHat);
	raytracingSphereShaderMaterial.uniforms.focusDistance.value = focusDistance;

	raytracingSphereShaderMaterial.uniforms.cosAlpha.value = Math.cos(alpha);
	raytracingSphereShaderMaterial.uniforms.sinAlpha.value = Math.sin(alpha);

}

/** create raytracing phere */
function addRaytracingSphere() {
	const videoFeedUTexture = new THREE.VideoTexture( videoFeedU );
	const videoFeedETexture = new THREE.VideoTexture( videoFeedE );
	videoFeedUTexture.colorSpace = THREE.SRGBColorSpace;
	videoFeedETexture.colorSpace = THREE.SRGBColorSpace;

	// create arrays of random numbers (as GLSL is rubbish at doing random numbers)
	let randomNumbersX = [];
	let randomNumbersY = [];
	// make the first random number 0 in both arrays, meaning the 0th ray starts from the centre of the aperture
	randomNumbersX.push(0);
	randomNumbersY.push(0);
	// fill in the rest of the array with random numbers
	let i=1;
	do {
		// create a new pairs or random numbers (x, y) such that x^2 + y^2 <= 1
		let x = 2*Math.random()-1;	// random number between -1 and 1
		let y = 2*Math.random()-1;	// random number between -1 and 1
		if(x*x + y*y <= 1) {
			// (x,y) lies within a circle of radius 1
			//  add a new point to the array of points on the aperture
			randomNumbersX.push(x);
			randomNumbersY.push(y);
			i++;
		}
	} while (i < 100);

	// the sphere surrouning the camera in all directions
	const geometry = 
		new THREE.SphereGeometry( raytracingSphereRadius );
	raytracingSphereShaderMaterial = new THREE.ShaderMaterial({
		side: THREE.DoubleSide,
		// wireframe: true,
		uniforms: {
			visible: { value: true },
			perfectRotatorVisible: { value: true },
			period: { value: 0.001 },
			cosAlpha: { value: 1.0 },
			sinAlpha: { value: 0.0	},
			stretchFactor: { value: 1.0 },
			additionalF: { value: 1e10 },	// additional focal length of lenslet array (an additional lens in the same plane)
			componentModelMatrix: { value: new THREE.Matrix4() },
			componentModelMatrixInverse: { value: new THREE.Matrix4() },
			// centreOfComponent: { value: new THREE.Vector3(0, 0, 0) },	// principal point of lenslet (0, 0)
			// comparisonModelMatrix: { value: new THREE.Matrix4() },
			// comparisonModelMatrixInverse: { value: new THREE.Matrix4() },
			// centreOfPerfectRotator: { value: new THREE.Vector3(0, 0, 0) },
			centreOfObjectPlane: { value: centreOfObjectPlane },
			designViewPosition: { value: designViewPosition },
			radius: { value: 0.05 },	// radius
			backgroundTexture: { value: backgroundTexture },
			showVideoFeed: { value: true },
			videoFeedUTexture: { value: videoFeedUTexture }, 
			videoFeedETexture: { value: videoFeedETexture },
			tanHalfFovHU: { value: 1.0 },
			tanHalfFovVU: { value: 1.0 },
			tanHalfFovHE: { value: 1.0 },
			tanHalfFovVE: { value: 1.0 },
			halfWidthU: { value: 1.0 },
			halfHeightU: { value: 1.0 },
			halfWidthE: { value: 1.0 },
			halfHeightE: { value: 1.0 },
			videoDistance: { value: 10.0 },	// distance of the image of the video feed from the origin
			focusDistance: { value: 10.0 },
			cameraXHat: { value: new THREE.Vector3(1, 0, 0) },
			cameraYHat: { value: new THREE.Vector3(0, 1, 0) },
			apertureRadius: { value: 0.0 },
			randomNumbersX: { value: randomNumbersX },
			randomNumbersY: { value: randomNumbersY },
			noOfRays: { value: 1 }
		},
		vertexShader: `
			varying vec3 intersectionPoint;

			void main()	{
				// projectionMatrix, modelViewMatrix, position -> passed in from Three.js
				intersectionPoint = (modelMatrix * vec4(position, 1.0)).xyz;	// position.xyz;
				
  				gl_Position = projectionMatrix
					* modelViewMatrix
					* vec4(position, 1.0);
			}
		`,
		fragmentShader: `
			precision highp float;

			#define PI 3.1415926538

			varying vec3 intersectionPoint;
			
			// the wedge array
			uniform mat4 componentModelMatrix;
			uniform mat4 componentModelMatrixInverse;
			uniform bool visible;
			uniform bool perfectRotatorVisible;
			uniform float cosAlpha;	// cos of rotation angle
			uniform float sinAlpha;	// sin of rotation angle
			uniform float stretchFactor;
			uniform float period;	// period of array
			uniform float additionalF;	// additional focal length (an additional lens in the same plane)
			// uniform vec3 centreOfComponent;	// centre of wedge array
			// uniform mat4 comparisonModelMatrix;
			// uniform mat4 comparisonModelMatrixInverse;
			// uniform vec3 centreOfPerfectRotator;
			uniform float radius;	// radius of wedge array
			uniform vec3 centreOfObjectPlane;
			uniform vec3 designViewPosition;

			// background
			uniform sampler2D backgroundTexture;
			uniform bool showVideoFeed;

			// video feed from user-facing camera
			uniform sampler2D videoFeedUTexture;
			uniform float halfWidthU;
			uniform float halfHeightU;

			// video feed from environment-facing camera
			uniform sampler2D videoFeedETexture;
			uniform float halfWidthE;
			uniform float halfHeightE;

			// the camera's wide aperture
			uniform float videoDistance;
			uniform float focusDistance;
			uniform int noOfRays;
			uniform vec3 cameraXHat;
			uniform vec3 cameraYHat;
			uniform float apertureRadius;
			uniform float randomNumbersX[100];
			uniform float randomNumbersY[100];

			// rotate the 2D vector v by the angle alpha (in radians)
			// from https://gist.github.com/yiwenl/3f804e80d0930e34a0b33359259b556c
			vec2 rotate(vec2 v, float cosAlpha, float sinAlpha) {
				mat2 m = mat2(cosAlpha, sinAlpha, -sinAlpha, cosAlpha);
				return m * v;
			}

			// returns a 3D matrix that corresponds to rotation by <angle> (in radians)
			// around <axis>
			// source: https://github.com/dmnsgn/glsl-rotate/blob/main/rotation-3d.glsl
			// use: rotation3d(axis, angle) * v gives the vector v, rotated by angle around axis
			mat4 rotation3d(vec3 axis, float angle) {
				axis = normalize(axis);
				float s = sin(angle);
				float c = cos(angle);
				float oc = 1.0 - c;

				return mat4(
					oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,  0.0,
					oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,  0.0,
					oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c,           0.0,
					0.0,                                0.0,                                0.0,                                1.0
				);
			}

			// propagate the ray starting at position p and with direction d to the plane z = z0, provided that plane
			// is in the ray's "forward" direction;
			// p becomes the point where the ray intersects p;
			// returns true or false depending on whether the intersection point is forwards or backwards along the ray
			bool propagateForwardToZPlane(
				inout vec3 p, 
				vec3 d, 
				float z0
			) {
				// calculate the z distance from the ray start position to the plane
				float deltaZ = z0 - p.z;

				// is the intersection with the plane in the ray's "forward" direction?
				if(d.z*deltaZ > 0.0) {
					// the intersection is in the forward direction; advance the ray to the plane
					p += d/d.z*deltaZ;	// set p to the intersection point with the plane
					return true;
				}
				return false;
			}

			// propagate the ray starting at position p and with direction d to the plane r.nHat = n0, provided that plane
			// is in the ray's "forward" direction;
			// p becomes the point where the ray intersects p;
			// returns true or false depending on whether the intersection point is forwards or backwards along the ray
			// nHat must be normalised
			bool propagateForwardToPlane(
				inout vec3 p,
				vec3 d,
				vec3 nHat,
				float n0
			) {
				// calculate the distance in the nHat direction from the ray start position to the plane
				float deltaN = n0 - dot(p, nHat);

				// is the intersection with the plane in the ray's "forward" direction?
				float dN = dot(d, nHat);
				if(dN*deltaN > 0.0) {
					// the intersection is in the forward direction; advance the ray to the plane
					p += d/dN*deltaN;	// set p to the intersection point with the plane
					return true;
				}
				return false;
			}

			// Calculate the light-ray direction after transmission through a lens or lens hologram.
			// d is the incident light-ray direction;
			// pixy is a 2D vector containing the transverse (x, y) components of the vector I-P,
			// i.e. the vector from the principal point P to the intersection point I;
			// f is the focal length;
			// returns the outgoing light-ray direction
			vec3 lensDeflect(vec3 d, vec2 pixy, float f) {
				bool idealLenses = true;
				if(idealLenses) {
					// ideal thin lens
					// "normalise" the direction such that the magnitude of the z component is 1
					vec3 d1 = d/abs(d.z);

					// the 3D deflected direction comprises the transverse components and a z component of magnitude 1
					// and the same sign as d.z
					return vec3(d1.xy - pixy/f, d1.z);
				} else {
					// lens hologram
					// normalise d
					vec3 dN = d/length(d);
					// transverse components of the outgoing light-ray direction
					vec2 dxy = dN.xy - pixy/f;
	
					// from the transverse direction, construct a 3D vector by setting the z component such that the length
					// of the vector is 1
					return vec3(dxy, sign(d.z)*sqrt(1.0 - dot(dxy, dxy)));
					}
			}

			// Pass the current ray (start point p, direction d, brightness factor b) through (or around) a lens.
			// The (ideal thin) lens, of focal length f, is in a z plane through centreOfLens.
			// It is circular, with the given radius, centred on centreOfLenss.
			void passThroughLens(
				inout vec3 p, 
				inout vec3 d, 
				inout vec4 b,
				vec3 centreOfLens, 
				float radius,
				float focalLength
			) {
				if(propagateForwardToZPlane(p, d, centreOfLens.z)) {
					// there is an intersection with the plane of this lens in the ray's forward direction

					// does the intersection point lie within the radius?
					vec2 pixy = p.xy - centreOfLens.xy;
					float r2 = dot(pixy, pixy);
					if(r2 < radius*radius) {
						// the intersection point lies inside the radius, so the lens does something to the ray

						// deflect the light-ray direction accordingly and make sure that the sign of the z component remains the same
						lensDeflect(d, pixy, focalLength);

						// lower the brightness factor, giving the light a blue tinge
						b *= vec4(0.9, 0.9, 0.99, 1);
					} 
				}
			}

			// find the u-coordinate of the centre of the nearest pixel in a 1D pixel array of period uPeriod
			float findPixelCentreCoordinate(float u, float uPeriod) {
				return uPeriod*floor(u/uPeriod+0.5);
			}

			// Find the centre of the nearest pixel in a rectangular pixel array centred on the origin.
			// r is a 2D vector containing the transverse components of the vector from the array centre to the
			// intersection point;
			// alpha is the angle (in radians) by which the array is rotated
			// period1 and period2 are the periods in the (rotated) x and y directions
			vec2 findPixelCentre(vec2 r, float cosAlpha, float sinAlpha, float period1, float period2) {
				vec2 rr = rotate(r, cosAlpha, -sinAlpha);
				vec2 lensletCentreST = vec2(
					findPixelCentreCoordinate(rr.x, period1),
					findPixelCentreCoordinate(rr.y, period2)
				);
				return rotate(lensletCentreST, cosAlpha, sinAlpha);
			}

			vec3 lensletArrayDeflect(vec3 d, vec3 intersectionPoint, float cosAlpha, float sinAlpha, float period, float focalLength) {
				vec2 p00ixy = intersectionPoint.xy;
				vec2 pxy = findPixelCentre(p00ixy, cosAlpha, sinAlpha, period, period);
				// light-ray direction after array
				return lensDeflect(d, intersectionPoint.xy - pxy, focalLength);
			}

			// Pass the current ray (start point p, direction d, brightness factor b) through (or around) a lenslet array.
			// The lenslet array is in a z plane through the origin; it is simulated as ideal thin lenses.
			// The lenslets, all of focal length f, are arranged in a square array of the given period, 
			// rotated by alpha around the z axis.  The component is circular, with the given radius, centred on the origin.
			void passThroughLensletArray(
				inout vec3 p, 
				inout vec3 d, 
				inout vec4 b,
				float radius,
				float cosAlpha,
				float sinAlpha,
				float period,
				float lensletsF,
				float overallF
			) {
				if(propagateForwardToZPlane(p, d, 0.)) {
					// there is an intersection with the plane of this array in the ray's forward direction
					
					// does the intersection point lie with in the radius?
					vec2 pixy = p.xy;
					float r2 = dot(pixy, pixy);	// length squared of vector r
					if(r2 < radius*radius) {
						// the intersection point lies inside the radius, so the component does something to the ray

						// deflect the light-ray direction accordingly 

						d = lensletArrayDeflect(d, p, cosAlpha, sinAlpha, period, lensletsF);
						d = lensDeflect(d, pixy, overallF);

						// lower the brightness factor, giving the light a blue tinge
						b *= vec4(0.9, 0.9, 0.99, 1);
					} 
				} // else b *= vec4(0.99, 0.9, 0.9, 1);	// this shouldn't happen -- give the light a red tinge
			}

			// deflection by an array of wedges in the plane z=0, centred on the origin
			vec3 zWedgeArrayDeflect(
				vec3 d, 
				vec3 intersectionPoint
			) {
				vec2 cixy = intersectionPoint.xy;	// - centreOfComponent.xy;

				// find the pixel centre
				vec3 p = vec3(findPixelCentre(cixy, 1., 0., period, period), 0);	// + centreOfComponent.xy, centreOfComponent.z);
				
				// calculate the incident and outgoing directions for which the pixel is designed
				vec3 i = normalize(p - designViewPosition);	// incident direction: from the design view position to the pixel centre
				// the x and y coordinates of the point where the undeviated ray from the designViewPosition through the pixel centre intersects the object plane
				vec2 undeviatedIntersectionPointXY = (p + centreOfObjectPlane.z*i/i.z).xy;	// (p + (centreOfObjectPlane.z - centreOfComponent.z)*i/i.z).xy;
				// the point where the deviated ray from the designViewPosition through the pixel centre intersects the object plane
				vec3 deviatedIntersectionPoint = 
					centreOfObjectPlane 
					+ vec3(rotate(undeviatedIntersectionPointXY-centreOfObjectPlane.xy, cosAlpha, -sinAlpha) / stretchFactor, 0.);
				// outgoing direction
				vec3 o = normalize(deviatedIntersectionPoint - p);

				// calculate the change deltaDxy in the normalised light-ray direction this pixel introduces
				vec2 deltaDXY = o.xy - i.xy;

				// calculate the transverse part of the outgoing light-ray direction
				vec2 t = normalize(d).xy + deltaDXY;

				// is the length of the transverse part > 1?
				float tSquared = dot(t, t);
				if(tSquared > 1.) {
					// yes, the outgoing ray is evanescent, so reflect it
					return vec3(d.xy, -d.z);
				}

				// complement the transverse part with the longitudinal component and return 
				return vec3(t, sign(d.z)*sqrt(1. - tSquared));
			}

			// Deflect a light ray transmitted through a view-rotating, pixellated, array of phase-holographic wedges.
			// This component is centred at the origin.	// <centreOfComponent>
			vec3 wedgeArrayDeflect(
				vec3 d, 
				vec3 intersectionPoint
			) {
				// the transverse directions are given by cameraXHat and cameraYHat,
				// the longitudinal direction is given by viewDirection

				vec3 ci = intersectionPoint;	// - centreOfComponent;
				vec2 cixy = vec2( dot(ci, cameraXHat), dot(ci, cameraYHat) );

				// find the pixel centre
				vec3 p = vec3(findPixelCentre(cixy, 1., 0., period, period),0);	// + centreOfComponent.xy, centreOfComponent.z);
				
				// calculate the incident and outgoing directions for which the pixel is designed
				vec3 i = normalize(p - designViewPosition);	// incident direction: from the design view position to the pixel centre
				// the x and y coordinates of the point where the undeviated ray from the designViewPosition through the pixel centre intersects the object plane
				vec2 undeviatedIntersectionPointXY = (p + centreOfObjectPlane.z*i/i.z).xy;	// (p + (centreOfObjectPlane.z - centreOfComponent.z)*i/i.z).xy;
				// the point where the deviated ray from the designViewPosition through the pixel centre intersects the object plane
				vec3 deviatedIntersectionPoint = 
					centreOfObjectPlane 
					+ vec3(rotate(undeviatedIntersectionPointXY-centreOfObjectPlane.xy, cosAlpha, -sinAlpha) / stretchFactor, 0.);
				// outgoing direction
				vec3 o = normalize(deviatedIntersectionPoint - p);

				// calculate the change deltaDxy in the normalised light-ray direction this pixel introduces
				vec2 deltaDXY = o.xy - i.xy;

				// calculate the transverse part of the outgoing light-ray direction
				vec2 t = normalize(d).xy + deltaDXY;

				// is the length of the transverse part > 1?
				float tSquared = dot(t, t);
				if(tSquared > 1.) {
					// yes, the outgoing ray is evanescent, so reflect it
					return vec3(d.xy, -d.z);
				}

				// complement the transverse part with the longitudinal component and return 
				return vec3(t, sign(d.z)*sqrt(1. - tSquared));
			}

			void perfectRotatorDeflect(
				inout vec3 p,
				inout vec3 d
			) {
				p = vec3(rotate(p.xy, cosAlpha, sinAlpha), 0.);
				d = vec3(rotate(d.xy, cosAlpha, sinAlpha), d.z);
			}

			// Pass the current ray (start point p, direction d, brightness factor b) through (or around) the wedge array.
			// The (holographic) wedge array is in the plane z=0
			// The wedges form a square array of the given period.
			// The component is circular, with the given radius, centred on the origin.
			// Each wedge is designed to make ...
			void passThroughWedgeArray(
				inout vec3 p, 
				inout vec3 d, 
				inout vec4 b,
				bool wedgeArrayVisible,
				bool perfectRotatorVisible
			) {
				p = (componentModelMatrixInverse*vec4(p,1)).xyz;
				d = (componentModelMatrixInverse*vec4(d,0)).xyz;
				if(propagateForwardToZPlane(p, d, 0.)) {	// centreOfComponent.z)) {
					// there is an intersection with the plane of this array in the ray's forward direction
					
					// does the intersection point lie within the radius?
					vec2 cixy = p.xy;	// - centreOfComponent.xy;
					float r2 = dot(cixy, cixy);	// length squared of vector r
					if(r2 < radius*radius) {
						// the intersection point lies inside the radius, so the component does something to the ray

						// deflect the light-ray direction accordingly 

						if(wedgeArrayVisible) {
							d = zWedgeArrayDeflect(d, p);

							// lower the brightness factor, giving the light a blue tinge
							b *= vec4(0.9, 0.9, 0.99, 1);
						}
						if(perfectRotatorVisible) {
							perfectRotatorDeflect(p, d);

							// lower the brightness factor, giving the light a blue tinge
							b *= vec4(0.9, 0.9, 0.99, 1);
						}
						// d = lensDeflect(d, cixy, overallF);
					} 
				} // else b *= vec4(0.99, 0.9, 0.9, 1);	// this shouldn't happen -- give the light a red tinge
				p = (componentModelMatrix*vec4(p,1)).xyz;
				d = (componentModelMatrix*vec4(d,0)).xyz;
			}

/* 			void passThroughPerfectRotator(
				inout vec3 p, 
				inout vec3 d, 
				inout vec4 b
			) {
				p = (comparisonModelMatrixInverse*vec4(p,1)).xyz;
				d = (comparisonModelMatrixInverse*vec4(d,0)).xyz;

				// if(propagateForwardToZPlane(p, d, 0.)) {
					// there is an intersection with the plane of this component in the ray's forward direction
					
					// does the intersection point lie within the radius?
					vec2 cixy = p.xy;	// - centreOfPerfectRotator.xy;
					float r2 = dot(cixy, cixy);	// length squared of vector r
					if(r2 < radius*radius) {
						// the intersection point lies inside the radius, so the component does something to the ray

						// change the light ray accordingly 
						perfectRotatorDeflect(p, d);

						// lower the brightness factor, giving the light a blue tinge
						b *= vec4(0.9, 0.9, 0.99, 1);
					} 
				// } // else b *= vec4(0.99, 0.9, 0.9, 1);	// this shouldn't happen -- give the light a red tinge

				p = (comparisonModelMatrix*vec4(p,1)).xyz;
				d = (comparisonModelMatrix*vec4(d,0)).xyz;
			}
 */
			// propagate the ray to the plane of the video feed, which is a z-distance <videoDistance> away,
			// and return either the color of the corresponding video-feed texel or the background color
			vec4 getColorOfVideoFeed(
				inout vec3 p, 
				vec3 d, 
				vec4 b,
				float videoFeedZ,
				sampler2D videoFeedTexture,
				float halfWidth,
				float halfHeight,
				vec4 backgroundColor
			) {
				// is the intersection in the ray's forward direction?
				if(propagateForwardToZPlane(p, d, videoFeedZ)) {
					// does the ray intersect the image?
					if((abs(p.x) < halfWidth) && (abs(p.y) < halfHeight))
						// yes, the ray intersects the image; take the pixel colour from the camera's video feed
						return texture2D(videoFeedTexture, vec2(0.5+0.5*p.x/halfWidth, 0.5+0.5*p.y/halfHeight));
					else 
						// the ray doesn't intersect the image
						return backgroundColor;
				}
			}

			vec4 getColorOfBackground(
				vec3 d
			) {
				float l = length(d);
				float phi = atan(d.z, d.x) + PI;
				float theta = acos(d.y/l);
				return texture2D(backgroundTexture, vec2(mod(phi/(2.*PI), 1.0), 1.-theta/PI));
			}

			void main() {
				// first calculate the point this pixel is focussed on, which is in a z plane a distance
				// <focusDistance> in  front of the camera, in the direction from the camera's aperture centre to the intersection point
				vec3 dF = intersectionPoint - cameraPosition;
				vec3 focusPosition = cameraPosition + focusDistance/abs(dF.z)*dF;

				// trace <noOfRays> rays
				gl_FragColor = vec4(0, 0, 0, 0);
				vec4 color;
				for(int i=0; i<noOfRays; i++) {
					// the current ray start position, a random point on the camera's circular aperture
					vec3 p = cameraPosition + apertureRadius*randomNumbersX[i]*cameraXHat + apertureRadius*randomNumbersY[i]*cameraYHat;
	
					// first calculate the current light-ray direction:
					// the ray first passes through focusPosition and then p,
					// so the "backwards" ray direction from the camera to the intersection point is
					//   d = focusPosition - p
					vec3 d = focusPosition - p;
					d = dF.z/d.z*d;
	
					// current brightness factor; this will multiply the colour at the end
					vec4 b = vec4(1.0, 1.0, 1.0, 1.0);
	
					if(d.z < 0.0) {
						// the ray is travelling "forwards", in the (-z) direction;
						// pass first through the array, then to environment-facing video feed
						passThroughWedgeArray(p, d, b, visible, perfectRotatorVisible);
						// passThroughPerfectRotator(p, d, b);
						if(showVideoFeed) color = getColorOfVideoFeed(p, d, b, -videoDistance, videoFeedETexture, halfWidthE, halfHeightE, vec4(1, 1, 1, 1.0));
						else color = getColorOfBackground(d);
					} else {
						// the ray is travelling "backwards", in the (+z) direction;
						// pass first through the array, then to user-facing video feed
						passThroughWedgeArray(p, d, b, visible, perfectRotatorVisible);
						// passThroughPerfectRotator(p, d, b);
						if(showVideoFeed) color = getColorOfVideoFeed(p, d, b, videoDistance, videoFeedUTexture, halfWidthU, halfHeightU, vec4(1, 1, 1, 1.0));
						else color = getColorOfBackground(d);
					}
		
					// finally, multiply by the brightness factor and add to gl_FragColor
					gl_FragColor += b*color;
				}
					
				gl_FragColor /= float(noOfRays);
			}
		`
	});
	raytracingSphere = new THREE.Mesh( geometry, raytracingSphereShaderMaterial ); 
	scene.add( raytracingSphere );
}

// see https://github.com/mrdoob/three.js/blob/master/examples/webgl_animation_skinning_additive_blending.html
function createGUI() {
	// const 
	gui = new GUI();
	// gui.hide();

	const params = {
		// 'Visible': raytracingSphereShaderMaterial.uniforms.visible.value,
		visible: function() {
			raytracingSphereShaderMaterial.uniforms.visible.value = !raytracingSphereShaderMaterial.uniforms.visible.value;
			visibleControl.name( visible2String() );
		},
		perfectRotatorVisible: function() {
			raytracingSphereShaderMaterial.uniforms.perfectRotatorVisible.value = !raytracingSphereShaderMaterial.uniforms.perfectRotatorVisible.value;
			perfectRotatorVisibleControl.name( perfectRotatorVisible2String() );
		},
		inFrontOfCamera: function() { 
			inFrontOfCamera = !inFrontOfCamera;
			inFrontOfCameraControl.name( inFrontOfCamera2String() );
		},
		ipd: 1000*ipd,
		'Component radius (cm)': 100.*raytracingSphereShaderMaterial.uniforms.radius.value,
		'tan<sup>-1</sup>(additional <i>F</i><sub>1</sub>)': Math.atan(raytracingSphereShaderMaterial.uniforms.additionalF.value),
		'Period, <i>p</i> (mm)': 1000.*raytracingSphereShaderMaterial.uniforms.period.value,
		'Rotation angle (&deg;)': alpha / Math.PI * 180.,
		'Stretch factor': raytracingSphereShaderMaterial.uniforms.stretchFactor.value,
		componentDistance: componentDistance,
		componentY: componentY,
		makeEyeLevel: function() { 
			componentY = camera.position.y;
			// componentYControl.setValue(componentY);
			GUIMesh.position.y = componentY - 1.5;
		},
		'Horiz. FOV (&deg;)': fovScreen,
		'Aperture radius (mm)': 1000.*raytracingSphereShaderMaterial.uniforms.apertureRadius.value,
		'tan<sup>-1</sup>(focus. dist.) (m)': Math.atan(focusDistance),
		'No of rays': noOfRays,
		'Env.-facing cam. (&deg;)': fovVideoFeedE,
		'User-facing cam. (&deg;)': fovVideoFeedU,
		'tan<sup>-1</sup>(video dist.)': Math.atan(raytracingSphereShaderMaterial.uniforms.videoDistance.value),
		'Point (virtual) cam. forward (in -<b>z</b> direction)': pointForward,
		'Show/hide info': toggleInfoVisibility,
		'Restart video streams': function() { 
			recreateVideoFeeds(); 
			postStatus("Restarting video stream");
		},
		background: function() {
			background = (background + 1) % 6;
			loadBackgroundImage();
			backgroundControl.name( background2String() );	
		},
		vrControlsVisible: function() {
			GUIMesh.visible = !GUIMesh.visible;
			vrControlsVisibleControl.name( guiMeshVisible2String() );
		},
	}

	// gui.add( params, 'Visible').onChange( (v) => { raytracingSphereShaderMaterial.uniforms.visible.value = v; } );
	visibleControl = gui.add( params, 'visible' ).name( visible2String() );
	perfectRotatorVisibleControl = gui.add( params, 'perfectRotatorVisible' ).name( perfectRotatorVisible2String() );
	inFrontOfCameraControl = gui.add( params, 'inFrontOfCamera' ).name( inFrontOfCamera2String() );
	gui.add( params, 'ipd', 0, 100, 1 ).onChange( (v) => { ipd = v/1000; } ).name( 'ipd (mm)' );
	gui.add( params, 'Component radius (cm)', 0, 10).onChange( (r) => { raytracingSphereShaderMaterial.uniforms.radius.value = r/100.; } );
	gui.add( params, 'Period, <i>p</i> (mm)', 0.1, 10).onChange( (p) => { raytracingSphereShaderMaterial.uniforms.period.value = p/1000.; } );
	// gui.add( params, 'tan<sup>-1</sup>(additional <i>F</i><sub>1</sub>)', -0.5*Math.PI, 0.5*Math.PI).onChange( (f) => { raytracingSphereShaderMaterial.uniforms.additionalF.value = Math.tan(f); } );
	gui.add( params, 'Rotation angle (&deg;)', -90, 90).onChange( (a) => { alpha = a/180.0*Math.PI; } );
	gui.add( params, 'Stretch factor', 0.1, 10 ).onChange( (m) => { raytracingSphereShaderMaterial.uniforms.stretchFactor.value = m; } );
	gui.add( params, 'componentDistance', 0.001, 0.2, 0.001 ).name( 'distance' ).onChange( (d) => { componentDistance = d; } );
	componentYControl = gui.add( params, 'componentY',  0, 3, 0.001).name( "<i>y</i><sub>centre</sub>" ).onChange( (y) => { componentY = y; } );
	gui.add( params, 'makeEyeLevel' ).name( 'Move to eye level' );
	
	gui.add( params, 'Point (virtual) cam. forward (in -<b>z</b> direction)');
	gui.add( params, 'Horiz. FOV (&deg;)', 10, 170, 1).onChange( setScreenFOV );
	gui.add( params, 'Aperture radius (mm)', 0, 10).onChange( (r) => { raytracingSphereShaderMaterial.uniforms.apertureRadius.value = r/1000.; } );
	gui.add( params, 'tan<sup>-1</sup>(focus. dist.) (m)', 
		//Math.atan(0.1), 
		-0.5*Math.PI,
		0.5*Math.PI
	).onChange( (a) => { focusDistance = Math.tan(a); } );
	gui.add( params, 'No of rays', 1, 100, 1).onChange( (n) => { noOfRays = n; } );

	backgroundControl = gui.add( params, 'background' ).name( background2String() );

	gui.add( params, 'Restart video streams');
	gui.add( params, 'Env.-facing cam. (&deg;)', 10, 170, 1).onChange( (fov) => { fovVideoFeedE = fov; updateUniforms(); });   
	gui.add( params, 'User-facing cam. (&deg;)', 10, 170, 1).onChange( (fov) => { fovVideoFeedU = fov; updateUniforms(); });  
	gui.add( params, 'tan<sup>-1</sup>(video dist.)', Math.atan(0.1), 0.5*Math.PI).onChange( (a) => { raytracingSphereShaderMaterial.uniforms.videoDistance.value = Math.tan(a); } );
	gui.add( params, 'Show/hide info');

	// if(renderer.xr.enabled) {
		vrControlsVisibleControl = gui.add( params, 'vrControlsVisible' );

		// create the GUI mesh at the end to make sure that it includes all controls
		GUIMesh = new HTMLMesh( gui.domElement );
/* 		GUIMesh.position.x = 0;
		GUIMesh.position.y = 1;	// componentY - 1.5;
		GUIMesh.position.z = -0.4;
		GUIMesh.rotation.x = 0;	// -Math.PI/4;
		GUIMesh.scale.setScalar( 2 );
		scene.add( GUIMesh );	
 */	
		GUIMesh.visible = false;
		vrControlsVisibleControl.name( guiMeshVisible2String() );	// this can be called only after GUIMesh has been created
	// }
}

function visible2String() {
	return 'Rotator '+(raytracingSphereShaderMaterial.uniforms.visible.value?'visible':'hidden');
}

function perfectRotatorVisible2String() {
	return 'Perfect rotator '+ (raytracingSphereShaderMaterial.uniforms.perfectRotatorVisible.value?'visible':'hidden');
}

function guiMeshVisible2String() {
	return 'VR controls '+(GUIMesh.visible?'visible':'hidden');
}

function inFrontOfCamera2String() {
	return 'Rotators ' + (inFrontOfCamera?'in front of Camera':'at origin');
}

function addXRInteractivity() {
	// see https://github.com/mrdoob/three.js/blob/master/examples/webxr_vr_sandbox.html

	// the two hand controllers

	const geometry = new THREE.BufferGeometry();
	geometry.setFromPoints( [ new THREE.Vector3( 0, 0, 0 ), new THREE.Vector3( 0, 0, - 5 ) ] );

	const controller1 = renderer.xr.getController( 0 );
	controller1.add( new THREE.Line( geometry ) );
	scene.add( controller1 );

	const controller2 = renderer.xr.getController( 1 );
	controller2.add( new THREE.Line( geometry ) );
	scene.add( controller2 );

	//

	const controllerModelFactory = new XRControllerModelFactory();

	const controllerGrip1 = renderer.xr.getControllerGrip( 0 );
	controllerGrip1.add( controllerModelFactory.createControllerModel( controllerGrip1 ) );
	scene.add( controllerGrip1 );

	const controllerGrip2 = renderer.xr.getControllerGrip( 1 );
	controllerGrip2.add( controllerModelFactory.createControllerModel( controllerGrip2 ) );
	scene.add( controllerGrip2 );

	//

	const group = new InteractiveGroup( renderer, camera );
	group.listenToPointerEvents( renderer, camera );
	group.listenToXRControllerEvents( controller1 );
	group.listenToXRControllerEvents( controller2 );
	scene.add( group );

	// place this below the resonator
	// GUIMesh = new HTMLMesh( gui.domElement );
 	GUIMesh.position.x = 0;
	GUIMesh.position.y = 1;	// componentY - 1.5;
	GUIMesh.position.z = -0.4;
	GUIMesh.rotation.x = 0;	// -Math.PI/4;
	GUIMesh.scale.setScalar( 2 );
	group.add( GUIMesh );	

 	// console.log('XR interactivity added');
}


function addDragControls() {
	let objects = [];
	objects.push(GUIMesh);

	dragControls = new DragControls( objects, camera, renderer.domElement );

	// add event listener to highlight dragged objects
	dragControls.addEventListener( 'dragstart', function ( event ) {
		event.object.material.emissive.set( 0xaaaaaa );
	} );

	dragControls.addEventListener( 'dragend', function ( event ) {
		event.object.material.emissive.set( 0x000000 );
	} );
}

function loadBackgroundImage() {
	// first free up resources
	if(backgroundTexture) backgroundTexture.dispose();
	if(background != 0) {
		const textureLoader = new THREE.TextureLoader();
		// textureLoader.crossOrigin = "Anonymous";

		let filename;
		switch (background) { 
			case 1:
				filename = '360-180 Glasgow University - Western Square.jpg';	// https://www.flickr.com/photos/pano_philou/1041580126
				break;
			case 2: 
				filename = '360-180 Glasgow University - Eastern Square.jpg';	// https://www.flickr.com/photos/pano_philou/1141564032
				break;
			case 3: 
				filename = 'Mugdock Woods 6 Milngavie Scotland Equirectangular.jpg';	// https://www.flickr.com/photos/gawthrop/3485817556
				break;
			case 4: 
				filename = 'Bluebells_13_Mugdock_Woods_Scotland-Equirectangular.jpg';	// https://www.flickr.com/photos/gawthrop/49889830418
				break;
			case 5: 
				filename = '360-180 The Glencoe Pass And The Three Sisters.jpg';	// https://www.flickr.com/photos/pano_philou/1140758031
				break;
			default:
				filename = '???';
				// 'Tower_University_Glasgow_Scotland-Equirectangular.jpg'	// https://www.flickr.com/photos/gawthrop/49890100126
				// 'Saddle_05_Arran_Scotland-Equirectangular.jpg'	// https://www.flickr.com/photos/gawthrop/49889356918
			}

		backgroundTexture = textureLoader.load(filename);
	}
}

function background2String() {
	switch (background) { 
	case 0: return 'Live video feed(s)';
	case 1: return 'Glasgow University, West Quadrangle';	// '360-180 Glasgow University - Western Square.jpg'	// https://www.flickr.com/photos/pano_philou/1041580126
	case 2: return 'Glasgow University, East Quadrangle';	// '360-180 Glasgow University - Eastern Square.jpg'	// https://www.flickr.com/photos/pano_philou/1141564032
	case 3: return 'Mugdock';	// 'Mugdock Woods 6 Milngavie Scotland Equirectangular.jpg'	// https://www.flickr.com/photos/gawthrop/3485817556
	case 4: return 'Mugdock bluebells';	// 'Bluebells_13_Mugdock_Woods_Scotland-Equirectangular.jpg'	// https://www.flickr.com/photos/gawthrop/49889830418
	case 5: return 'Glencoe';	// '360-180 The Glencoe Pass And The Three Sisters.jpg'	// https://www.flickr.com/photos/pano_philou/1140758031
	default: return 'Undefined';		
		// 'Tower_University_Glasgow_Scotland-Equirectangular.jpg'	// https://www.flickr.com/photos/gawthrop/49890100126
		// 'Saddle_05_Arran_Scotland-Equirectangular.jpg'	// https://www.flickr.com/photos/gawthrop/49889356918
	}
}

function createVideoFeeds() {
	// create the video stream for the user-facing camera first, as some devices (such as my iPad), which have both cameras,
	// but can (for whatever reason) only have a video feed from one at a time, seem to go with the video stream that was
	// created last, and as the standard view is looking "forward" it is preferable to see the environment-facing camera.
	videoFeedU = document.getElementById( 'videoFeedU' );

	// see https://github.com/mrdoob/three.js/blob/master/examples/webgl_materials_video_webcam.html
	if ( navigator.mediaDevices && navigator.mediaDevices.getUserMedia ) {
		// user-facing camera
		const constraintsU = { video: { 
			// 'deviceId': cameraId,	// this could be the device ID selected 
			width: {ideal: 1280},	// {ideal: 10000}, 
			// height: {ideal: 10000}, 
			facingMode: {ideal: 'user'}
			// aspectRatio: { exact: width / height }
		} };
		navigator.mediaDevices.getUserMedia( constraintsU ).then( function ( stream ) {
			// apply the stream to the video element used in the texture
			videoFeedU.srcObject = stream;
			videoFeedU.play();

			videoFeedU.addEventListener("playing", () => {
				aspectRatioVideoFeedU = videoFeedU.videoWidth / videoFeedU.videoHeight;
				updateUniforms();
				postStatus(`User-facing(?) camera resolution ${videoFeedU.videoWidth} &times; ${videoFeedU.videoHeight}`);
			});
		} ).catch( function ( error ) {
			postStatus(`Unable to access user-facing camera/webcam (Error: ${error})`);
		} );
	} else {
		postStatus( 'MediaDevices interface, which is required for video streams from device cameras, not available.' );
	}

	videoFeedE = document.getElementById( 'videoFeedE' );

	// see https://github.com/mrdoob/three.js/blob/master/examples/webgl_materials_video_webcam.html
	if ( navigator.mediaDevices && navigator.mediaDevices.getUserMedia ) {
		// environment-facing camera
		const constraintsE = { video: { 
			// 'deviceId': cameraId,	// this could be the device ID selected 
			width: {ideal: 1280},	// {ideal: 10000}, 
			// height: {ideal: 10000}, 
			facingMode: {ideal: 'environment'}
			// aspectRatio: { exact: width / height }
		} };
		navigator.mediaDevices.getUserMedia( constraintsE ).then( function ( stream ) {
			// apply the stream to the video element used in the texture
			videoFeedE.srcObject = stream;
			videoFeedE.play();

			videoFeedE.addEventListener("playing", () => {
				aspectRatioVideoFeedE = videoFeedE.videoWidth / videoFeedE.videoHeight;
				updateUniforms();
				postStatus(`Environment-facing(?) camera resolution ${videoFeedE.videoWidth} &times; ${videoFeedE.videoHeight}`);
			});
		} ).catch( function ( error ) {
			postStatus(`Unable to access environment-facing camera/webcam (Error: ${error})`);
		} );
	} else {
		postStatus( 'MediaDevices interface, which is required for video streams from device cameras, not available.' );
	}
}

function addEventListenersEtc() {
	// handle device orientation
	// window.addEventListener("deviceorientation", handleOrientation, true);

	// handle window resize
	window.addEventListener("resize", onWindowResize, false);

	// handle screen-orientation (landscape/portrait) change
	screen.orientation.addEventListener( "change", recreateVideoFeeds );

	// share button functionality
	document.getElementById('takePhotoButton').addEventListener('click', takePhoto);

	// toggle fullscreen button functionality
	document.getElementById('fullscreenButton').addEventListener('click', toggleFullscreen);

	// info button functionality
	document.getElementById('infoButton').addEventListener('click', toggleInfoVisibility);

	// back button functionality
	document.getElementById('backButton').addEventListener('click', showLivePhoto);
	document.getElementById('backButton').style.visibility = "hidden";

	// share button
	document.getElementById('shareButton').addEventListener('click', share);
	document.getElementById('shareButton').style.visibility = "hidden";
	if(!(navigator.share)) document.getElementById('shareButton').src="./shareButtonUnavailable.png";
	// if(!(navigator.share)) document.getElementById('shareButton').style.opacity = 0.3;

	// delete button
	document.getElementById('deleteButton').addEventListener('click', deleteStoredPhoto);
	document.getElementById('deleteButton').style.visibility = "hidden";

	// hide the thumbnail for the moment
	document.getElementById('storedPhotoThumbnail').addEventListener('click', showStoredPhoto);
	document.getElementById('storedPhotoThumbnail').style.visibility = "hidden";
	document.getElementById('storedPhoto').addEventListener('click', showLivePhoto);
	document.getElementById('storedPhoto').style.visibility = "hidden";
	// showingStoredPhoto = false;
}

/**
 * @param {*} fov	The larger of the camera's horizontal and vertical FOV, in degrees
 * 
 * Set the larger FOV of the screen/window to fov.
 * 
 * Depending on the screen/window's FOV, fov is either the horizontal fov (if screen width > screen height)
 * or the vertical fov (if screen width < screen height).
 */
function setScreenFOV(fov) {
	fovScreen = fov;

	screenChanged();
}

// function swapArrays() {
// 	const visible3 = lensletArrayShaderMaterial.uniforms.visible1.value;
// 	const focalLength3 = lensletArrayShaderMaterial.uniforms.lensletsF1.value;
// 	const period3 = lensletArrayShaderMaterial.uniforms.period1.value;
// 	const alpha3 = lensletArrayShaderMaterial.uniforms.alpha1.value;

// 	lensletArrayShaderMaterial.uniforms.visible1.value = lensletArrayShaderMaterial.uniforms.visible2.value;
// 	lensletArrayShaderMaterial.uniforms.lensletsF1.value = lensletArrayShaderMaterial.uniforms.lensletsF2.value;
// 	lensletArrayShaderMaterial.uniforms.period1.value = lensletArrayShaderMaterial.uniforms.period2.value;
// 	lensletArrayShaderMaterial.uniforms.alpha1.value = lensletArrayShaderMaterial.uniforms.alpha2.value;

// 	lensletArrayShaderMaterial.uniforms.visible2.value = visible3;
// 	lensletArrayShaderMaterial.uniforms.lensletsF2.value = focalLength3;
// 	lensletArrayShaderMaterial.uniforms.period2.value = period3;
// 	lensletArrayShaderMaterial.uniforms.alpha2.value = alpha3;
// }

/** 
 * Reset the aspect ratio and FOV of the virtual cameras.
 * 
 * Call if the window size has changed (which also happens when the screen orientation changes)
 * or if camera's FOV has changed
 */
function screenChanged() {
	// alert(`new window size ${window.innerWidth} x ${window.innerHeight}`);

	// in case the screen size has changed
	if(renderer) renderer.setSize(window.innerWidth, window.innerHeight);

	// if the screen orientation changes, width and height swap places, so the aspect ratio changes
	let windowAspectRatio = window.innerWidth / window.innerHeight;
	camera.aspect = windowAspectRatio;

	// fovS is the screen's horizontal or vertical FOV, whichever is greater;
	// re-calculate the camera FOV, which is the *vertical* fov
	let verticalFOV;
	if(windowAspectRatio > 1.0) {
		// fovS is horizontal FOV; convert to get correct vertical FOV
		verticalFOV = 2.0*Math.atan(Math.tan(0.5*fovScreen*Math.PI/180.0)/windowAspectRatio)*180.0/Math.PI;
	} else {
		// fovS is already vertical FOV
		verticalFOV = fovScreen;
	}
	camera.fov = verticalFOV;

	// make sure the camera changes take effect
	camera.updateProjectionMatrix();
}

function  pointForward() {
	let r = camera.position.length();
	camera.position.x = 0;
	camera.position.y = 0;
	camera.position.z = r;
	orbitControls.update();
	postStatus('Pointing camera forwards (in -<b>z</b> direction)');
}

function onWindowResize() {
	screenChanged();
	postStatus(`window size ${window.innerWidth} &times; ${window.innerHeight}`);	// debug
}

// // see https://developer.mozilla.org/en-US/docs/Web/API/ScreenOrientation/change_event
function recreateVideoFeeds() {
	// stop current video streams...
	videoFeedE.srcObject.getTracks().forEach(function(track) { track.stop(); });
	videoFeedU.srcObject.getTracks().forEach(function(track) { track.stop(); });

	// ... and re-create new ones, hopefully of the appropriate size
	createVideoFeeds();
}

function addOrbitControls() {
	// controls

	orbitControls = new OrbitControls( camera, renderer.domElement );
	// controls = new OrbitControls( cameraOutside, renderer.domElement );
	orbitControls.listenToKeyEvents( window ); // optional

	//controls.addEventListener( 'change', render ); // call this only in static scenes (i.e., if there is no animation loop)
	orbitControls.addEventListener( 'change', cameraPositionChanged );

	orbitControls.enableDamping = false; // an animation loop is required when either damping or auto-rotation are enabled
	orbitControls.dampingFactor = 0.05;

	orbitControls.enablePan = true;
	orbitControls.enableZoom = true;

	orbitControls.maxPolarAngle = Math.PI;
}

function cameraPositionChanged() {
	postStatus(`Camera position (${camera.position.x.toPrecision(2)}, ${camera.position.y.toPrecision(2)}, ${camera.position.z.toPrecision(2)})`);
	// counter = 0;
	// keep the raytracing sphere centred on the camera position
	// raytracingSphere.position.copy(camera.position.clone());	// TODO this doesn't seem to work as intended!?
}

async function toggleFullscreen() {
	if (!document.fullscreenElement) {
		document.documentElement.requestFullscreen().catch((err) => {
			postStatus(
				`Error attempting to enable fullscreen mode: ${err.message} (${err.name})`,
			);
		});
		// allow screen orientation changes
		// screen.orientation.unlock();
	} else {
		document.exitFullscreen();
	}
}

function showStoredPhoto() {
	gui.hide();
	renderer.domElement.style.visibility = "hidden";
	document.getElementById('takePhotoButton').style.visibility = "hidden";
	// document.getElementById('changePositionButton').style.visibility = "hidden";
	document.getElementById('storedPhotoThumbnail').style.visibility = "hidden";
	document.getElementById('backButton').style.visibility = "visible";
	document.getElementById('shareButton').style.visibility = "visible";
	document.getElementById('deleteButton').style.visibility = "visible";
	document.getElementById('storedPhoto').style.visibility = "visible";
	showingStoredPhoto = true;

	postStatus('Showing stored photo, '+storedPhotoDescription);
}

function showLivePhoto() {
	gui.show();
	renderer.domElement.style.visibility = "visible";
	document.getElementById('takePhotoButton').style.visibility = "visible";
	// document.getElementById('changePositionButton').style.visibility = "visible";
	if(storedPhoto) document.getElementById('storedPhotoThumbnail').style.visibility = "visible";
	document.getElementById('backButton').style.visibility = "hidden";
	document.getElementById('shareButton').style.visibility = "hidden";
	document.getElementById('deleteButton').style.visibility = "hidden";
	document.getElementById('storedPhoto').style.visibility = "hidden";
	showingStoredPhoto = false;

	postStatus('Showing live image');
}

function deleteStoredPhoto() {
	bin.play();

	storedPhoto = null;

	showLivePhoto();

	postStatus('Stored photo deleted; showing live image');
}

function takePhoto() {
	try {
		click.play();

		storedPhoto = renderer.domElement.toDataURL('image/png');
		storedPhotoInfoString = getInfoString();

		storedPhotoDescription = name;
		// 
		document.getElementById('storedPhoto').src=storedPhoto;
		document.getElementById('storedPhotoThumbnail').src=storedPhoto;
		document.getElementById('storedPhotoThumbnail').style.visibility = "visible";
	
		postStatus('Photo taken; click thumbnail to view and share');
	} catch (error) {
		console.error('Error:', error);
	}	
}

async function share() {
	try {
		fetch(storedPhoto)
		.then(response => response.blob())
		.then(blob => {
			const file = new File([blob], name+storedPhotoDescription+'.png', { type: blob.type });

			// Use the Web Share API to share the screenshot
			if (navigator.share) {
				navigator.share({
					title: storedPhotoDescription,
					files: [file],
				});
			} else {
				postStatus('Sharing is not supported by this browser.');
			}	
		})
		.catch(error => {
			console.error('Error:', error);
			postStatus(`Error: ${error}`);
		});
	} catch (error) {
		console.error('Error:', error);
	}
}

/** 
 * Add a text field to the bottom left corner of the screen
 */
function createStatus() {
	// see https://stackoverflow.com/questions/15248872/dynamically-create-2d-text-in-three-js
	status.style.position = 'absolute';
	status.style.backgroundColor = "rgba(0, 0, 0, 0.3)";	// semi-transparent black
	status.style.color = "White";
	status.style.fontFamily = "Arial";
	status.style.fontSize = "9pt";
	postStatus("Welcome!");
	status.style.bottom = 0 + 'px';
	status.style.left = 0 + 'px';
	status.style.zIndex = 1;
	document.body.appendChild(status);	
}

function postStatus(text) {
	status.innerHTML = '&nbsp;'+text;
	// console.log('status: '+text);

	// show the text only for 3 seconds
	statusTime = new Date().getTime();
	setTimeout( () => { if(new Date().getTime() - statusTime > 2999) status.innerHTML = '&nbsp;'+name+', University of Glasgow, <a href="https://github.com/jkcuk/'+name+'">https://github.com/jkcuk/'+name+'</a>' }, 3000);
}

function getInfoString() {
	return `Array<br>` +
		`&nbsp;&nbsp;Visible `+ (raytracingSphereShaderMaterial.uniforms.visible.value?'&check;':'&cross;')+`<br>` +
		`&nbsp;&nbsp;Period = ${raytracingSphereShaderMaterial.uniforms.period.value.toPrecision(4)}<br>` +
		`&nbsp;&nbsp;Rotation angle = ${(alpha*180.0/Math.PI).toPrecision(4)}&deg;<br>` +
		`&nbsp;&nbsp;Radius = ${raytracingSphereShaderMaterial.uniforms.radius.value.toPrecision(4)}<br>` +
		// `&nbsp;&nbsp;Centre of wedge array = (${raytracingSphereShaderMaterial.uniforms.centreOfComponent.value.x.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.centreOfComponent.value.y.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.centreOfComponent.value.z.toPrecision(4)})<br>` +
		`&nbsp;&nbsp;Focal length of additional lens in same plane = ${raytracingSphereShaderMaterial.uniforms.additionalF.value.toPrecision(4)}<br>` +		
		`Video feeds<br>` +
		`&nbsp;&nbsp;Distance from origin = ${raytracingSphereShaderMaterial.uniforms.videoDistance.value.toPrecision(4)}<br>` +	// (user-facing) camera
		`&nbsp;&nbsp;Horizontal fields of view (when seen from the origin)<br>` +
		`&nbsp;&nbsp;&nbsp;&nbsp;User-facing camera = ${fovVideoFeedU.toPrecision(4)}&deg;<br>` +	// (user-facing) camera
		`&nbsp;&nbsp;&nbsp;&nbsp;Environment-facing camera = ${fovVideoFeedE.toPrecision(4)}&deg;<br>` +	// (environment-facing) camera
		`Virtual camera<br>` +
		`&nbsp;&nbsp;Position = (${camera.position.x.toPrecision(4)}, ${camera.position.y.toPrecision(4)}, ${camera.position.z.toPrecision(4)})<br>` +
		`&nbsp;&nbsp;Horiz. FOV = ${fovScreen.toPrecision(4)}<br>` +
		`&nbsp;&nbsp;Aperture radius = ${raytracingSphereShaderMaterial.uniforms.apertureRadius.value.toPrecision(4)}<br>` +
		`&nbsp;&nbsp;Focussing distance = ${focusDistance.toPrecision(4)}<br>` +
		`&nbsp;&nbsp;Number of rays = ${noOfRays}`
		// `cameraXHat = (${raytracingSphereShaderMaterial.uniforms.cameraXHat.value.x.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.cameraXHat.value.y.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.cameraXHat.value.z.toPrecision(4)})<br>` +
		// `cameraYHat = (${raytracingSphereShaderMaterial.uniforms.cameraYHat.value.x.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.cameraYHat.value.y.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.cameraYHat.value.z.toPrecision(4)})`
		;
		console.log("*");
}

function refreshInfo() {
	if(showingStoredPhoto) setInfo( storedPhotoInfoString );
	else setInfo( getInfoString() );

	if(info.style.visibility == "visible") setTimeout( refreshInfo , 100);	// refresh again a while
}

/** 
 * Add a text field to the top left corner of the screen
 */
function createInfo() {
	// see https://stackoverflow.com/questions/15248872/dynamically-create-2d-text-in-three-js
	info.style.position = 'absolute';
	info.style.backgroundColor = "rgba(0, 0, 0, 0.3)";	// semi-transparent black
	info.style.color = "White";
	info.style.fontFamily = "Arial";
	info.style.fontSize = "9pt";
	info.innerHTML = "-- nothing to show (yet) --";
	info.style.top = 60 + 'px';
	info.style.left = 0 + 'px';
	info.style.zIndex = 1;
	document.body.appendChild(info);
	info.style.visibility = "hidden";
}

function setInfo(text) {
	info.innerHTML = text;
	console.log('info: '+text);
	// // show the text only for 3 seconds
	// infoTime = new Date().getTime();
	// setTimeout( () => { if(new Date().getTime() - infoTime > 2999) info.innerHTML = `` }, 3000);
	// info.style.visibility = "visible";
}

function toggleInfoVisibility() {
	switch(info.style.visibility) {
		case "visible":
			info.style.visibility = "hidden";
			break;
		case "hidden":
		default:
			info.style.visibility = "visible";
			refreshInfo();
	}
}