use std::{cmp::Ordering, fmt::Debug};
use num_integer::Integer;

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

impl crate::SquareRoot for Num {
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

    fn from_str_radix(_str: &str, _radix: u32) -> Result<Self, Self::FromStrRadixErr> {
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

    fn rem(self, _rhs: Self) -> Self::Output {
        todo!()
    }
}