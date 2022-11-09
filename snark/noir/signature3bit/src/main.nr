use dep::std;

fn main(
    priv_key : [Field; 3]
    ) {
    // We make it first 
    let mut accum = [0,0,0,0];
    for i in 0..3 {
        accum[i+1] = priv_key[i] + 2 * accum[i];
    };
    for i in 0..4 {
        if i > 0 {
            constrain priv_key[i-1] * priv_key[i-1] == priv_key[i-1];
        }
    };
}