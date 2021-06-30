extern crate num_traits;

use std::{cmp::Ordering, fmt::Debug};
use euclid::*;
use bevy::prelude::*;
use euclid::Vector3D;
use num_integer::Integer;
use num_traits::Signed;

#[derive(Copy, Clone, Debug)]
struct Num {
    // a*sqrt(b)/c
    a: i128,
    b: u128,
    c: u128,
}

impl Num {
    pub fn new(a: i128, b: u128, c: u128) -> Self {
        Num{a,b,c}.normalize()
    }

    pub fn normalize(mut self) -> Self {
        self.reduce_sqrt();
        self.reduce_fraction();

        self
    }

    pub fn reduce_sqrt(&mut self) {
        let mut i = 2u128;
        while i*i<=self.b {
            if self.b%(i*i)==0 {
                self.b/=i*i;
                self.a*=i as i128;
            }else{
                i+=1;
            }
        }
    }

    pub fn reduce_fraction(&mut self) {
        let divisor = if self.a.is_negative() {
            (self.a.abs() as u128).gcd(&self.c)
        }else{
            (self.a as u128).gcd(&self.c)
        };
        self.a/=divisor as i128;
        self.c/=divisor;
    }
}

impl SquareRoot for Num {
    fn sqrt(self) -> Self {
        let mut x = self;
        assert!(x.a>0);
        x.reduce_sqrt();
        x.reduce_fraction();
        assert!(x.b==1);
        Num::new(1,x.a as u128*x.c,x.c)
    }
}

impl Into<f64> for Num {
    fn into(self) -> f64 {
        self.a as f64 * (self.b as f64).sqrt() / self.c as f64
    }
}

impl From<f32> for Num {
    fn from(a: f32) -> Self {
        assert!(a.is_finite());
        let two_to_23 = 8_388_608;
        let x = (a*two_to_23 as f32) as i128;

        Num::new(x,1,two_to_23 as u128)
    }
}

impl PartialOrd for Num {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.a.is_negative() ^ other.a.is_negative() {
            //one is negative
            if self.a.is_negative() {
                Some(Ordering::Less)
            }else{
                Some(Ordering::Greater)
            }
        }else if self.a.is_negative() {
            //both are negative -> flip the comparison
            let left = (self.a*self.a) as u128*other.c*other.c*self.b;
            let right = (other.a*other.a) as u128*self.c*self.c*other.b;
            right.partial_cmp(&left)
        }else{
            //both are positive
            let left = (self.a*self.a) as u128*other.c*other.c*self.b;
            let right = (other.a*other.a) as u128*self.c*self.c*other.b;
            left.partial_cmp(&right)
        }
    }
}

impl PartialEq for Num {
    fn eq(&self, other: &Self) -> bool {
        let x = self.normalize();
        let y = other.normalize();
        x.a==y.a && x.b==y.b && x.c==y.c
    }
}

impl num_traits::Num for Num {
    type FromStrRadixErr=();

    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        todo!()
    }
}

impl num_traits::Zero for Num {
    fn zero() -> Self {
        Num::new(0, 1, 1)
    }

    fn is_zero(&self) -> bool {
        self.a==0 || self.b==0
    }
}

impl num_traits::One for Num {
    fn one() -> Self {
        Num::new(1, 1, 1)
    }

    fn is_one(&self) -> bool {
        let x = self.normalize();
        x.a==1 && x.b==1 && x.c==1
    }
}

impl std::ops::Add for Num {
    type Output=Num;

    fn add(self, rhs: Self) -> Self::Output {
        if self.b==rhs.b {
            // a*sqrt(b)/c + d*sqrt(b)/e = sqrt(b)*(a/c + d/e) = sqrt(b)*(ae+dc)/ce
            let mut num = Num {
                a: self.a*rhs.c as i128 + rhs.a*self.c as i128,
                b: self.b,
                c: self.c*self.c,
            };
            num.reduce_fraction();
            num
        }else{
            todo!()
        }
    }
}

impl std::ops::Sub for Num {
    type Output=Num;

    fn sub(self, rhs: Self) -> Self::Output {
        if self.b==rhs.b {
            // a*sqrt(b)/c + d*sqrt(b)/e = sqrt(b)*(a/c + d/e) = sqrt(b)*(ae+dc)/ce
            let mut num = Num {
                a: self.a*rhs.c as i128 - rhs.a*self.c as i128,
                b: self.b,
                c: self.c*self.c,
            };
            num.reduce_fraction();
            num
        }else{
            todo!()
        }
    }
}

impl std::ops::Mul for Num {
    type Output=Num;

    fn mul(self, rhs: Self) -> Self::Output {
        Num {
            a: self.a*rhs.a,
            b: self.b*rhs.b,
            c: self.c*rhs.c,
        }.normalize()
    }
}

impl std::ops::Div for Num {
    type Output=Num;

    fn div(self, rhs: Self) -> Self::Output {
        if self.a.is_negative() ^ rhs.a.is_negative() {
            // result negative
            Num {
                a: -1*self.a.abs()*rhs.c as i128,
                b: self.b*rhs.b,
                c: self.c*rhs.a.abs() as u128*rhs.b,
            }.normalize()
        }else{
            //result positive
            Num {
                a: self.a.abs()*rhs.c as i128,
                b: self.b*rhs.b,
                c: self.c*rhs.a.abs() as u128*rhs.b,
            }.normalize()
        }
    }
}

impl std::ops::Rem for Num {
    type Output=Num;

    fn rem(self, rhs: Self) -> Self::Output {
        todo!()
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
        let c = difference.square_length()-self.radius_sq;
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
        let vector = object.center-self.center;
        let center_distance_squared = vector.square_length();
        let center_distance = center_distance_squared.sqrt();
        let vector_normalized = vector/center_distance;
        /*
        r1^2-x^2=r^2
        r2^2-(d-x)^2=r^2

        r1^2-x^2=r2^2-(d-x)^2
        r1^2-x^2=r2^2-d^2+2dx-x^2
        d^2+r1^2-r2^2=2dx
        (d^2+r1^2-r2^2)/2d=x
        */

        let circle_center_distance_from_self_center = (center_distance_squared+self.radius_sq-object.radius_sq)/(T::from(2.0)*center_distance);
        let circle_center = self.center+vector_normalized*circle_center_distance_from_self_center;
        let circle_radius_squared = self.radius_sq-circle_center_distance_from_self_center*circle_center_distance_from_self_center;
        Box::new([Circle{
            center: circle_center,
            radius_sq: circle_radius_squared,
            normal: vector,     // or normalized vector?
        }])
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
    radius_sq: T,
}

pub struct Circle<T> where T: Copy + num_traits::Num {
    center: Point3D<T,UnknownUnit>,
    radius_sq: T,
    normal: Vector3D<T,UnknownUnit>,    // or normalized vector?
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

    //test Num
    let x = Num::from(6.125);
    dbg!(x);

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
        radius_sq:25.0,
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

    //drawing ball
    commands.spawn_bundle(PbrBundle {
        mesh: meshes.add(Mesh::from(shape::Icosphere { radius: ball.radius_sq.sqrt() as f32, subdivisions: 5 })),
        material: materials.add(Color::rgb(0.8, 0.7, 0.6).into()),
        transform: Transform::from_translation(ball.center.to_tuple().into()),
        ..Default::default()
    });

    let ball2: Ball<f32> = Ball{
        center:(4.0,2.0,0.0).into(),
        radius_sq:16.0,
    };

    let circle = ball.intersections(&ball2);

    if circle.len()>0 {
        //drawing circle
        let circle = &circle[0];
        commands.spawn_bundle(PbrBundle {
            mesh: meshes.add(Mesh::from(shape::Torus { radius: circle.radius_sq.sqrt(), ring_radius: 0.1, subdivisions_segments: 32, subdivisions_sides: 24 })),
            material: materials.add(Color::BLUE.into()),
            transform: Transform {translation: circle.center.to_tuple().into(), rotation: Quat::from_rotation_arc(Vec3::Y,circle.normal.normalize().to_tuple().into()), scale: Vec3::ONE},
            ..Default::default()
        });

    }

    //drawing second ball
    commands.spawn_bundle(PbrBundle {
        mesh: meshes.add(Mesh::from(shape::Icosphere { radius: ball2.radius_sq.sqrt() as f32, subdivisions: 5 })),
        material: materials.add(Color::rgb(0.8, 0.7, 0.6).into()),
        transform: Transform::from_translation(ball2.center.to_tuple().into()),
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
