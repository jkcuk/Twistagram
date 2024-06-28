
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

let cameraLensDistance = 0.02;
let alpha = 0.0;

// camera with wide aperture
let focusDistance = 1e8;
let noOfRays = 1;

let raytracingSphereRadius = 10.0;	// what this code does should be independent of this value, but isn't (?)

// "internal" variables
let raytracingSphere;
let raytracingSphereShaderMaterial;

let scene;
let renderer;
let videoFeedU, videoFeedE;	// feeds from user/environment-facing cameras
let camera;
let controls;

// the menu
let gui;
let GUIMesh;
let vrControlsVisibleControl;

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
	camera = new THREE.PerspectiveCamera( fovScreen, windowAspectRatio, 0.01*raytracingSphereRadius, 2*raytracingSphereRadius );
	camera.position.z = cameraLensDistance;
	screenChanged();
	
	renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
	renderer.setPixelRatio(window.devicePixelRatio);
	renderer.setSize( window.innerWidth, window.innerHeight );
	document.body.appendChild( renderer.domElement );
	// document.getElementById('livePhoto').appendChild( renderer.domElement );

	createVideoFeeds();

	addRaytracingSphere();

	// user interface

	addEventListenersEtc();

	addOrbitControls();

	// the controls menu
	createGUI();

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
	requestAnimationFrame( animate );

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

	// create the points on the aperture

	// create basis vectors for the camera's clear aperture
	let viewDirection = new THREE.Vector3();
	let apertureBasisVector1 = new THREE.Vector3();
	let apertureBasisVector2 = new THREE.Vector3();
	camera.getWorldDirection(viewDirection);
	// if(counter < 10) console.log(`viewDirection = (${viewDirection.x.toPrecision(2)}, ${viewDirection.y.toPrecision(2)}, ${viewDirection.z.toPrecision(2)})`);

	if((viewDirection.x == 0.0) && (viewDirection.y == 0.0)) {
		// viewDirection is along z direction
		apertureBasisVector1.crossVectors(viewDirection, new THREE.Vector3(1, 0, 0)).normalize();
	} else {
		// viewDirection is not along z direction
		apertureBasisVector1.crossVectors(viewDirection, new THREE.Vector3(0, 0, 1)).normalize();
	}
	apertureBasisVector2.crossVectors(viewDirection, apertureBasisVector1).normalize();

	raytracingSphereShaderMaterial.uniforms.noOfRays.value = noOfRays;
	raytracingSphereShaderMaterial.uniforms.apertureXHat.value.copy(apertureBasisVector1);
	raytracingSphereShaderMaterial.uniforms.apertureYHat.value.copy(apertureBasisVector2);
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
			period: { value: 0.001 },
			cosAlpha: { value: 1.0 },
			sinAlpha: { value: 0.0	},
			stretchFactor: { value: 1.0 },
			additionalF: { value: 1e10 },	// additional focal length of lenslet array (an additional lens in the same plane)
			centreOfWedgeArray: { value: new THREE.Vector3(0, 0, 0) },	// principal point of lenslet (0, 0)
			centreOfObjectPlane: { value: new THREE.Vector3(0, 0, -10) },
			designViewPosition: { value: new THREE.Vector3(0, 0, 0.02) },
			radius: { value: 0.05 },	// radius
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
			apertureXHat: { value: new THREE.Vector3(1, 0, 0) },
			apertureYHat: { value: new THREE.Vector3(0, 1, 0) },
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

			varying vec3 intersectionPoint;
			
			// the wedge array
			uniform bool visible;
			uniform float cosAlpha;	// cos of rotation angle
			uniform float sinAlpha;	// sin of rotation angle
			uniform float stretchFactor;
			uniform float period;	// period of array
			uniform float additionalF;	// additional focal length (an additional lens in the same plane)
			uniform vec3 centreOfWedgeArray;	// centre of wedge array
			uniform float radius;	// radius of wedge array
			uniform vec3 centreOfObjectPlane;
			uniform vec3 designViewPosition;

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
			uniform vec3 apertureXHat;
			uniform vec3 apertureYHat;
			uniform float apertureRadius;
			uniform float randomNumbersX[100];
			uniform float randomNumbersY[100];

			// rotate the 2D vector v by the angle alpha (in radians)
			// from https://gist.github.com/yiwenl/3f804e80d0930e34a0b33359259b556c
			vec2 rotate(vec2 v, float cosAlpha, float sinAlpha) {
				mat2 m = mat2(cosAlpha, sinAlpha, -sinAlpha, cosAlpha);
				return m * v;
			}

			// propagate the ray starting at position p and with direction d to the plane z = z0, providing that plane
			// is in the ray's "forward" direction;
			// p becomes the point where the ray intersects p;
			// returns true or false depending on whether the intersection point is forwards or backwards along the ray
			bool propagateForwardToZPlane(
				inout vec3 p, 
				vec3 d, 
				float z0
			) {
				// calculate the z distance from the ray start position to the array
				float deltaZ = z0 - p.z;

				// is the intersection with the plane in the ray's "forward" direction?
				if(d.z*deltaZ > 0.0) {
					// the intersection is in the forward direction; advance the ray to the plane
					p += d/d.z*deltaZ;	// set p to the intersection point with the plane
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

			// Find the centre of the nearest pixel in a rectangular pixel array.
			// r is a 2D vector containing the transverse components of the vector from pp00 to the
			// intersection point;
			// pp00 is the centre of the pixel with indices (0, 0) -- sort of the array centre;
			// alpha is the angle (in radians) by which the array is rotated
			// period1 and period2 are the periods in the (rotated) x and y directions
			vec2 findPixelCentre(vec2 r, vec2 pp00, float cosAlpha, float sinAlpha, float period1, float period2) {
				vec2 rr = rotate(r, cosAlpha, -sinAlpha);
				vec2 lensletCentreST = vec2(
					findPixelCentreCoordinate(rr.x, period1),
					findPixelCentreCoordinate(rr.y, period2)
				);
				return rotate(lensletCentreST, cosAlpha, sinAlpha) + pp00;
			}

			vec3 lensletArrayDeflect(vec3 d, vec3 intersectionPoint, vec3 principalPoint00, float cosAlpha, float sinAlpha, float period, float focalLength) {
				vec2 p00ixy = intersectionPoint.xy - principalPoint00.xy;
				vec2 pxy = findPixelCentre(p00ixy, principalPoint00.xy, cosAlpha, sinAlpha, period, period);
				// light-ray direction after array
				return lensDeflect(d, intersectionPoint.xy - pxy, focalLength);
			}

			// Pass the current ray (start point p, direction d, brightness factor b) through (or around) a lenslet array.
			// The lenslet array is in a z plane through centreOfWedgeArray; it is simulated as ideal thin lenses.
			// The lenslets, all of focal length f, are arranged in a square array of the given period, 
			// rotated by alpha around the z axis.  The component is circular, with the given radius, centred on centreOfWedgeArray.
			void passThroughLensletArray(
				inout vec3 p, 
				inout vec3 d, 
				inout vec4 b,
				vec3 centreOfArray, 
				float radius,
				float cosAlpha,
				float sinAlpha,
				float period,
				float lensletsF,
				float overallF
			) {
				if(propagateForwardToZPlane(p, d, centreOfArray.z)) {
					// there is an intersection with the plane of this array in the ray's forward direction
					
					// does the intersection point lie with in the radius?
					vec2 pixy = p.xy - centreOfArray.xy;
					float r2 = dot(pixy, pixy);	// length squared of vector r
					if(r2 < radius*radius) {
						// the intersection point lies inside the radius, so the component does something to the ray

						// deflect the light-ray direction accordingly 

						d = lensletArrayDeflect(d, p, centreOfArray, cosAlpha, sinAlpha, period, lensletsF);
						d = lensDeflect(d, pixy, overallF);

						// lower the brightness factor, giving the light a blue tinge
						b *= vec4(0.9, 0.9, 0.99, 1);
					} 
				} else b *= vec4(0.99, 0.9, 0.9, 1);	// this shouldn't happen -- give the light a red tinge
			}

			vec3 wedgeArrayDeflect(
				vec3 d, 
				vec3 intersectionPoint
				//vec3 centreOfWedgeArray, 
				//float cosAlpha, 
				//float sinAlpha, 
				//float stretchFactor,
				//vec3 centreOfObjectPlane,
				//vec3 designViewPosition,
				//float period
			) {
				vec2 cixy = intersectionPoint.xy - centreOfWedgeArray.xy;

				// find the pixel centre
				vec3 p = centreOfWedgeArray + vec3(findPixelCentre(cixy, centreOfWedgeArray.xy, 1., 0., period, period), 0.);
				
				// calculate the incident and outgoing directions for which the pixel is designed
				vec3 i = normalize(p - designViewPosition);	// incident direction: from the design view position to the pixel centre
				// the x and y coordinates of the point where the undeviated ray from the designViewPosition through the pixel centre intersects the object plane
				vec2 undeviatedIntersectionPointXY = (p + (centreOfObjectPlane.z - centreOfWedgeArray.z)*i/i.z).xy;
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

			// Pass the current ray (start point p, direction d, brightness factor b) through (or around) the wedge array.
			// The (holographic) wedge array is in a z plane through centreOfWedgeArray.
			// The wedges form a square array of the given period.
			// The component is circular, with the given radius, centred on centreOfWedgeArray.
			// Each wedge is designed to make ...
			void passThroughWedgeArray(
				inout vec3 p, 
				inout vec3 d, 
				inout vec4 b
				//vec3 centreOfWedgeArray, 
				//float radius,
				//float cosAlpha,
				//float sinAlpha,
				//float period
			) {
				if(propagateForwardToZPlane(p, d, centreOfWedgeArray.z)) {
					// there is an intersection with the plane of this array in the ray's forward direction
					
					// does the intersection point lie within the radius?
					vec2 cixy = p.xy - centreOfWedgeArray.xy;
					float r2 = dot(cixy, cixy);	// length squared of vector r
					if(r2 < radius*radius) {
						// the intersection point lies inside the radius, so the component does something to the ray

						// deflect the light-ray direction accordingly 

						d = wedgeArrayDeflect(d, p);
						// d = lensDeflect(d, cixy, overallF);

						// lower the brightness factor, giving the light a blue tinge
						b *= vec4(0.9, 0.9, 0.99, 1);
					} 
				} else b *= vec4(0.99, 0.9, 0.9, 1);	// this shouldn't happen -- give the light a red tinge
			}

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
					vec3 p = cameraPosition + apertureRadius*randomNumbersX[i]*apertureXHat + apertureRadius*randomNumbersY[i]*apertureYHat;
	
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
						if(visible) passThroughWedgeArray(p, d, b);
						color = getColorOfVideoFeed(p, d, b, -videoDistance, videoFeedETexture, halfWidthE, halfHeightE, vec4(1, 1, 1, 1.0));	// white
					} else {
						// the ray is travelling "backwards", in the (+z) direction;
						// pass first through the array, then to user-facing video feed
						if(visible) passThroughWedgeArray(p, d, b);
						color = getColorOfVideoFeed(p, d, b, videoDistance, videoFeedUTexture, halfWidthU, halfHeightU, vec4(1, 0, 0, 1.0));	// white
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
		'Visible': raytracingSphereShaderMaterial.uniforms.visible.value,
		'Component radius (cm)': 100.*raytracingSphereShaderMaterial.uniforms.radius.value,
		'tan<sup>-1</sup>(additional <i>F</i><sub>1</sub>)': Math.atan(raytracingSphereShaderMaterial.uniforms.additionalF.value),
		'Period, <i>p</i> (mm)': 1000.*raytracingSphereShaderMaterial.uniforms.period.value,
		'Rotation angle (&deg;)': alpha / Math.PI * 180.,
		'Stretch factor': raytracingSphereShaderMaterial.uniforms.stretchFactor.value,
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
		vrControlsVisible: function() {
			GUIMesh.visible = !GUIMesh.visible;
			vrControlsVisibleControl.name( guiMeshVisible2String() );
		},
	}

	gui.add( params, 'Visible').onChange( (v) => { raytracingSphereShaderMaterial.uniforms.visible.value = v; } );
	gui.add( params, 'Component radius (cm)', 0, 10).onChange( (r) => { raytracingSphereShaderMaterial.uniforms.radius.value = r/100.; } );
	gui.add( params, 'Period, <i>p</i> (mm)', 0.1, 10).onChange( (p) => { raytracingSphereShaderMaterial.uniforms.period.value = p/1000.; } );
	// gui.add( params, 'tan<sup>-1</sup>(additional <i>F</i><sub>1</sub>)', -0.5*Math.PI, 0.5*Math.PI).onChange( (f) => { raytracingSphereShaderMaterial.uniforms.additionalF.value = Math.tan(f); } );
	gui.add( params, 'Rotation angle (&deg;)', -90, 90).onChange( (a) => { alpha = a/180.0*Math.PI; } );
	gui.add( params, 'Stretch factor', 0.1, 10 ).onChange( (m) => { raytracingSphereShaderMaterial.uniforms.stretchFactor.value = m; } );
	gui.add( params, 'Point (virtual) cam. forward (in -<b>z</b> direction)');
	gui.add( params, 'Horiz. FOV (&deg;)', 10, 170, 1).onChange( setScreenFOV );
	gui.add( params, 'Aperture radius (mm)', 0, 10).onChange( (r) => { raytracingSphereShaderMaterial.uniforms.apertureRadius.value = r/1000.; } );
	gui.add( params, 'tan<sup>-1</sup>(focus. dist.) (m)', 
		//Math.atan(0.1), 
		-0.5*Math.PI,
		0.5*Math.PI
	).onChange( (a) => { focusDistance = Math.tan(a); } );
	gui.add( params, 'No of rays', 1, 100, 1).onChange( (n) => { noOfRays = n; } );

	gui.add( params, 'Restart video streams');
	gui.add( params, 'Env.-facing cam. (&deg;)', 10, 170, 1).onChange( (fov) => { fovVideoFeedE = fov; updateUniforms(); });   
	gui.add( params, 'User-facing cam. (&deg;)', 10, 170, 1).onChange( (fov) => { fovVideoFeedU = fov; updateUniforms(); });  
	gui.add( params, 'tan<sup>-1</sup>(video dist.)', Math.atan(0.1), 0.5*Math.PI).onChange( (a) => { raytracingSphereShaderMaterial.uniforms.videoDistance.value = Math.tan(a); } );
	gui.add( params, 'Show/hide info');

	if(renderer.xr.enabled) {
		vrControlsVisibleControl = gui.add( GUIParams, 'vrControlsVisible' );

		// create the GUI mesh at the end to make sure that it includes all controls
		GUIMesh = new HTMLMesh( gui.domElement );
		GUIMesh.visible = false;
		vrControlsVisibleControl.name( guiMeshVisible2String() );	// this can be called only after GUIMesh has been created
	}
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
	GUIMesh.position.y = resonatorY - 1.5;
	GUIMesh.position.z = -0.4;
	GUIMesh.rotation.x = -Math.PI/4;
	GUIMesh.scale.setScalar( 2 );
	group.add( GUIMesh );	
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
	controls.update();
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

	controls = new OrbitControls( camera, renderer.domElement );
	// controls = new OrbitControls( cameraOutside, renderer.domElement );
	controls.listenToKeyEvents( window ); // optional

	//controls.addEventListener( 'change', render ); // call this only in static scenes (i.e., if there is no animation loop)
	controls.addEventListener( 'change', cameraPositionChanged );

	controls.enableDamping = false; // an animation loop is required when either damping or auto-rotation are enabled
	controls.dampingFactor = 0.05;

	controls.enablePan = true;
	controls.enableZoom = true;

	controls.maxPolarAngle = Math.PI;
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
		`&nbsp;&nbsp;Centre of wedge array = (${raytracingSphereShaderMaterial.uniforms.centreOfWedgeArray.value.x.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.centreOfWedgeArray.value.y.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.centreOfWedgeArray.value.z.toPrecision(4)})<br>` +
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
		// `apertureXHat = (${raytracingSphereShaderMaterial.uniforms.apertureXHat.value.x.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.apertureXHat.value.y.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.apertureXHat.value.z.toPrecision(4)})<br>` +
		// `apertureYHat = (${raytracingSphereShaderMaterial.uniforms.apertureYHat.value.x.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.apertureYHat.value.y.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.apertureYHat.value.z.toPrecision(4)})`
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