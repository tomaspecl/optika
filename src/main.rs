extern crate num_traits;

use std::fmt::Debug;
use euclid::*;
use bevy::prelude::*;
use euclid::Vector3D;

struct Num {
    // a*sqrt(b)/c
    a: i128,
    b: u128,
    c: u128,
}

impl Num {
    pub fn new(a: i128, b: u128, c: u128) -> Self {
        Num{a,b,c}
    }
}

pub trait SquareRoot {
    fn sqrt(self) -> Self;
}

impl SquareRoot for f64 {
    fn sqrt(self) -> f64 {
        self.sqrt()
    }
}

impl SquareRoot for f32 {
    fn sqrt(self) -> f32 {
        self.sqrt()
    }
}

pub trait Parametric<T> where T: Copy + num_traits::Num {
    type Parameter;
    fn point(&self, parameter: Self::Parameter) -> Point3D<T,UnknownUnit>;
}

pub trait Intersect<O,T> where T: Copy + num_traits::Num {
    type Intersection;
    fn intersections(&self, object: &O) -> Box<[Self::Intersection]>;
}

pub trait Surface<T> where T: Copy + num_traits::Num {
    fn normal(&self, point: Point3D<T,UnknownUnit>) -> Vector3D<T,UnknownUnit>;
}

pub trait ParametricSurface<T>: Parametric<T> where T: Copy + num_traits::Num {
    fn normal(&self, parameter: Self::Parameter) -> Vector3D<T,UnknownUnit>;
}

pub trait Reflect<T>: Intersect<Line<T>,T> where T: Copy + num_traits::Num {
    fn reflect(&self, line: &Line<T>) -> Option<Line<T>>;
}

impl<T> Parametric<T> for Line<T> where T: Copy + num_traits::Num {
    type Parameter=T;

    fn point(&self, parameter: Self::Parameter) -> Point3D<T,UnknownUnit> {
        self.start+self.direction*parameter
    }
}

impl<T> Surface<T> for Ball<T> where T: Copy + num_traits::Num {
    fn normal(&self, point: Point3D<T,UnknownUnit>) -> Vector3D<T,UnknownUnit> {
        point-self.center
    }
}

impl<T> Intersect<Line<T>,T> for Ball<T> where T: Copy + num_traits::Num + SquareRoot + From<f32> + PartialOrd {
    type Intersection=Point3D<T,UnknownUnit>;

    // should return the point which is closer to the start of the line and then the farther one, also the parameter for the returned point has to be within the parameter bounds
    fn intersections(&self, object: &Line<T>) -> Box<[Self::Intersection]> {
        // (lx+t*sx-x0)^2 + (ly+t*sy-y0)^2 + (lz+t*sz-z0)^2 = r^2
        // (lx-x0)^2 + 2*(lx-x0)*t*sx + t^2*sx^2  +  (ly-y0)^2 + 2*(ly-y0)*t*sy + t^2*sy^2  +  (lz-z0)^2 + 2*(lz-z0)*t*sz + t^2*sz^2 = r^2
        // (lx-x0)^2 + (ly-y0)^2 + (lz-z0)^2  +  (2*sx*(lx-x0) + 2*sy*(ly-y0) + 2*sz*(lz-z0))*t  +  (sx^2 + sy^2 + sz^2)*t^2 = r^2
        // a = sx^2 + sy^2 + sz^2
        // b = 2*(sx*(lx-x0) + sy*(ly-y0) + sz*(lz-z0))
        // c = (lx-x0)^2 + (ly-y0)^2 + (lz-z0)^2 - r^2
        // D = b^2-4*a*c
        // t12 = (-b+-sqrt(D))/2a
        let difference = object.start-self.center;
        let a = object.direction.square_length();
        let b = T::from(2.0)*(object.direction.dot(difference));
        let c = difference.square_length()-self.radius*self.radius;
        let discriminant = b*b-T::from(4.0)*a*c;

        if discriminant>T::zero() {
            // t1 should be closer to 0 than t2
            let t1 = (T::from(-1.0)*b-discriminant.sqrt())/(T::from(2.0)*a);
            let t2 = (T::from(-1.0)*b+discriminant.sqrt())/(T::from(2.0)*a);
            if t1>=object.parameter_bounds.0 && t1 <= object.parameter_bounds.1 {
                let p1 = object.point(t1);
                if t2>=object.parameter_bounds.0 && t2 <= object.parameter_bounds.1 {
                    let p2 = object.point(t2);
                    Box::new([p1,p2])
                }else{
                    Box::new([p1])
                }
            }else{
                if t2>=object.parameter_bounds.0 && t2 <= object.parameter_bounds.1 {
                    let p2 = object.point(t2);
                    Box::new([p2])
                }else{
                    Box::new([])
                }
            }
        }else if discriminant==T::zero() {
            let t = (T::from(-1.0)*b)/(T::from(2.0)*a);
            if t>=T::zero() {
                let p = object.point(t);
                Box::new([p])
            }else{
                Box::new([])
            }
        }else{
            Box::new([])
        }
    }
}

impl<T> Intersect<Ball<T>,T> for Ball<T> where T: Copy + num_traits::Num + SquareRoot + From<f32> + PartialOrd {
    type Intersection=Circle<T>;

    fn intersections(&self, object: &Ball<T>) -> Box<[Self::Intersection]> {
        let center_distance = (self.center-object.center).square_length();
        todo!()
    }
}

#[derive(Debug)]
pub struct Line<T> where T: Copy + num_traits::Num {
    start: Point3D<T,UnknownUnit>,
    direction: Vector3D<T,UnknownUnit>,
    parameter_bounds: (T,T),
}

// (x-x0)^2 + (y-y0)^2 + (z-z0)^2 = r^2
#[derive(Debug)]
pub struct Ball<T> where T: Copy + num_traits::Num {
    center: Point3D<T,UnknownUnit>,
    radius: T,
}

pub struct Circle<T> where T: Copy + num_traits::Num {
    center: Point3D<T,UnknownUnit>,
    radius: T,
    normal: Vector3D<T,UnknownUnit>,
}

impl<T> Reflect<T> for Ball<T> where T: Copy + num_traits::Num + SquareRoot + From<f32> + PartialOrd {
    
    fn reflect(&self, line: &Line<T>) -> Option<Line<T>> {
        let point = self.intersections(line);
        if point.len()>0 {
            let point=point[0];
            let normal = self.normal(point);
            let reflected = line.direction-normal*T::from(2.0)*line.direction.dot(normal)/normal.square_length();
            Some(Line{
                start: point,
                direction: reflected,
                parameter_bounds: (T::zero(),T::one()/T::zero()),
            })
        }else{
            None
        }
        
    }
}

fn main() {
    App::build()
        .insert_resource(Msaa { samples: 4 })
        .add_plugins(DefaultPlugins)
        .add_startup_system(setup.system())
        .run();
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    /*// plane
    commands.spawn_bundle(PbrBundle {
        mesh: meshes.add(Mesh::from(shape::Plane { size: 5.0 })),
        material: materials.add(Color::rgb(0.3, 0.5, 0.3).into()),
        ..Default::default()
    });*/

    let line: Line<f32> = Line{
        direction:(1.0,2.0,-2.0).into(),
        start:(-5.0,-5.0,8.0).into(),
        parameter_bounds: (0.0,f32::INFINITY),
    };

    let ball: Ball<f32> = Ball{
        center:(0.0,0.0,0.0).into(),
        radius:5.0,
    };

    let now = std::time::Instant::now();
    let intersections = ball.intersections(&line);
    println!("elapsed: {}",now.elapsed().as_nanos());

    dbg!(&intersections);

    let now = std::time::Instant::now();
    let reflect = ball.reflect(&line);
    println!("elapsed: {}",now.elapsed().as_nanos());

    dbg!(&reflect);

    //drawing line
    commands.spawn_bundle(PbrBundle {
        mesh: meshes.add(Mesh::from(shape::Box{min_x:-0.05, max_x:0.05, min_y:-0.05, max_y:0.05, min_z:0.0, max_z:50.0})),
        material: materials.add(Color::rgb(0.8, 0.7, 0.6).into()),
        transform: Transform::from_translation(line.start.to_tuple().into()).looking_at((line.start-line.direction).to_tuple().into(), Vec3::Y),
        ..Default::default()
    });

    if let Some(line) = reflect {
        //drawing reflected line
        commands.spawn_bundle(PbrBundle {
            mesh: meshes.add(Mesh::from(shape::Box{min_x:-0.05, max_x:0.05, min_y:-0.05, max_y:0.05, min_z:0.0, max_z:50.0})),
            material: materials.add(Color::rgb(0.8, 0.7, 0.6).into()),
            transform: Transform::from_translation(line.start.to_tuple().into()).looking_at((line.start-line.direction).to_tuple().into(), Vec3::Y),
            ..Default::default()
        });
    }

    commands.spawn_bundle(PbrBundle {
        mesh: meshes.add(Mesh::from(shape::Icosphere { radius: ball.radius as f32, subdivisions: 5 })),
        material: materials.add(Color::rgb(0.8, 0.7, 0.6).into()),
        transform: Transform::from_translation(ball.center.to_tuple().into()),
        ..Default::default()
    });

    // light
    commands.spawn_bundle(LightBundle {
        transform: Transform::from_xyz(4.0, 8.0, 4.0),
        ..Default::default()
    });

    // camera
    commands.spawn_bundle(PerspectiveCameraBundle {
        transform: Transform::from_xyz(-2.0*5.0, 2.5*5.0, 5.0*5.0).looking_at(Vec3::ZERO, Vec3::Y),
        ..Default::default()
    });
}
