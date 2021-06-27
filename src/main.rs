use std::fmt::Debug;
use std::ops::*;
use std::vec;
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

pub trait Intersect<T> where T: Copy + num_traits::Num {
    type Object: Object<T>;
    fn intersections(self, object: Self::Object) -> Box<[Point3D<T,UnknownUnit>]>;
}

pub trait Object<T> where T: Copy + Add<T,Output=T> + Sub<T,Output=T> {
    type Parameter;
    fn point(self, pos: Self::Parameter) -> Point3D<T,UnknownUnit>;
    fn normal(self, pos: Self::Parameter) -> Vector3D<T,UnknownUnit>;
    fn intersections(self, line: &Line<T>) -> Box<[Self::Parameter]>;
}

impl<T: Copy + Add<T,Output=T> + Mul<T,Output=T> > Object<T> for Ball<T> {
    type Parameter=T;

    fn point(self, pos: Self::Parameter) -> Point3D<T,UnknownUnit> {

    }

    fn normal(self, pos: Self::Parameter) -> Vector3D<T,UnknownUnit> {
        self.point(pos)-self.center
    }

    /// Returns first and second intersection
    fn intersections(self, line: &Line<T>) -> Box<[Self::Parameter]> {
        // (lx+t*sx-x0)^2 + (ly+t*sy-y0)^2 + (lz+t*sz-z0)^2 = r^2
        // (lx-x0)^2 + 2*(lx-x0)*t*sx + t^2*sx^2  +  (ly-y0)^2 + 2*(ly-y0)*t*sy + t^2*sy^2  +  (lz-z0)^2 + 2*(lz-z0)*t*sz + t^2*sz^2 = r^2
        // (lx-x0)^2 + (ly-y0)^2 + (lz-z0)^2  +  (2*sx*(lx-x0) + 2*sy*(ly-y0) + 2*sz*(lz-z0))*t  +  (sx^2 + sy^2 + sz^2)*t^2 = r^2
        // a = sx^2 + sy^2 + sz^2
        // b = 2*(sx*(lx-x0) + sy*(ly-y0) + sz*(lz-z0))
        // c = (lx-x0)^2 + (ly-y0)^2 + (lz-z0)^2 - r^2
        // D = b^2-4*a*c
        // t12 = (-b+-sqrt(D))/2a
        let difference = line.start-self.center;
        let a = line.direction.square_length();
        let b = T::from(2.0)*(line.direction.dot(difference));
        let c = difference.square_length()-self.radius*self.radius;
        let discriminant = b*b-T::from(4.0)*a*c;

        if discriminant>=T::from(0.0) {
            let t1 = (-b+discriminant.sqrt())/(T::from(2.0)*a);
            let t2 = (-b-discriminant.sqrt())/(T::from(2.0)*a);

            let p1 = line.point(t1);
            let p2 = line.point(t2);

            // TODO: is the order correct?
            Box::new([t2,t1])
        }else if discriminant==T::from(0.0) {
            let t = -b/(T::from(2.0)*a);
            Box::new([t])
        }else{
            Box::new([])
        }
    }
}

#[derive(Debug)]
pub struct Line<T> where T: Copy + Add<T,Output=T> {
    start: Point3D<T,UnknownUnit>,
    direction: Vector3D<T,UnknownUnit>,
}

impl<T> Line<T> where T: Copy + Mul<T,Output=T> + Add<T,Output=T> {
    pub fn point(&self, parameter: T) -> Point3D<T,UnknownUnit> {
        self.start+self.direction*parameter
    }
}

// (x-x0)^2 + (y-y0)^2 + (z-z0)^2 = r^2
#[derive(Debug)]
struct Ball<T> where T: Copy + Mul + Add {
    center: Point3D<T,UnknownUnit>,
    radius: T,
}

impl<T> Ball<T> where T: Copy + Debug + From<f32> + Mul<T,Output=T> + Div<T,Output=T> + Add<T,Output=T> + Sub<T,Output=T> + Neg<Output=T> + SquareRoot + std::cmp::PartialOrd {
    
    pub fn reflect(&self, line: &Line<T>) -> Option<Line<T>> {
        let point = self.intersection(line)?.0;
        let normal = point-self.center;
        let reflected = line.direction-normal*T::from(2.0)*line.direction.dot(normal)/normal.square_length();
        Some(Line{
            start: point,
            direction: reflected,
        })
    }
    
    pub fn intersection(&self, line: &Line<T>) -> Option<(Point3D<T,UnknownUnit>, Point3D<T,UnknownUnit>)>{
        // (lx+t*sx-x0)^2 + (ly+t*sy-y0)^2 + (lz+t*sz-z0)^2 = r^2
        // (lx-x0)^2 + 2*(lx-x0)*t*sx + t^2*sx^2  +  (ly-y0)^2 + 2*(ly-y0)*t*sy + t^2*sy^2  +  (lz-z0)^2 + 2*(lz-z0)*t*sz + t^2*sz^2 = r^2
        // (lx-x0)^2 + (ly-y0)^2 + (lz-z0)^2  +  (2*sx*(lx-x0) + 2*sy*(ly-y0) + 2*sz*(lz-z0))*t  +  (sx^2 + sy^2 + sz^2)*t^2 = r^2
        // a = sx^2 + sy^2 + sz^2
        // b = 2*(sx*(lx-x0) + sy*(ly-y0) + sz*(lz-z0))
        // c = (lx-x0)^2 + (ly-y0)^2 + (lz-z0)^2 - r^2
        // D = b^2-4*a*c
        // t12 = (-b+-sqrt(D))/2/a
        let difference = line.start-self.center;
        let a = line.direction.square_length();
        let b = T::from(2.0)*(line.direction.dot(difference));
        let c = difference.square_length()-self.radius*self.radius;
        let discriminant = b*b-T::from(4.0)*a*c;

        if discriminant>=T::from(0.0) {
            let t1 = (-b+discriminant.sqrt())/(T::from(2.0)*a);
            let t2 = (-b-discriminant.sqrt())/(T::from(2.0)*a);

            let p1 = line.point(t1);
            let p2 = line.point(t2);

            // TODO: is the order correct?
            Some((p2,p1))
        }else{
            None
        }
        
    }
}

fn main() {
    let line = Line{
        direction:(1.0,0.0,0.0).into(),
        start:(-5.0,0.0,-5.0).into(),
    };

    let ball = Ball{
        center:(0.0,0.0,0.0).into(),
        radius:5.0,
    };
    let now = std::time::Instant::now();
    let intersections = ball.intersection(&line);
    println!("elapsed: {}",now.elapsed().as_nanos());

    dbg!(intersections);

    let now = std::time::Instant::now();
    let reflect = ball.reflect(&line);
    println!("elapsed: {}",now.elapsed().as_nanos());

    dbg!(reflect);

    println!("equal: {}", intersections.unwrap().0==intersections.unwrap().1);

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
    };
    let ball: Ball<f32> = Ball{
        center:(0.0,0.0,0.0).into(),
        radius:5.0,
    };

    commands.spawn_bundle(PbrBundle {
        mesh: meshes.add(Mesh::from(shape::Icosphere { radius: ball.radius as f32, subdivisions: 5 })),
        material: materials.add(Color::rgb(0.8, 0.7, 0.6).into()),
        transform: Transform::from_translation(ball.center.to_tuple().into()),
        ..Default::default()
    });

    
    //line
    commands.spawn_bundle(PbrBundle {
        mesh: meshes.add(Mesh::from(shape::Box{min_x:-0.05, max_x:0.05, min_y:-0.05, max_y:0.05, min_z:0.0, max_z:50.0})),
        material: materials.add(Color::rgb(0.8, 0.7, 0.6).into()),
        transform: Transform::from_translation(line.start.to_tuple().into()).looking_at((line.start-line.direction).to_tuple().into(), Vec3::Y),
        ..Default::default()
    });

    let line = ball.reflect(&line).unwrap();

    //reflected line
    commands.spawn_bundle(PbrBundle {
        mesh: meshes.add(Mesh::from(shape::Box{min_x:-0.05, max_x:0.05, min_y:-0.05, max_y:0.05, min_z:0.0, max_z:50.0})),
        material: materials.add(Color::rgb(0.8, 0.7, 0.6).into()),
        transform: Transform::from_translation(line.start.to_tuple().into()).looking_at((line.start-line.direction).to_tuple().into(), Vec3::Y),
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
