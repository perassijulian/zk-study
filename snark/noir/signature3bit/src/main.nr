use dep::std;

fn main(
    priv_key : [Field; 3],
    ) {
    for i in 0..3 {
        constrain priv_key[i] * priv_key[i] == priv_key[i];
    };
    let mut accum = [0,0,0,0];
    for i in 0..3 {
        accum[i+1] = priv_key[i] + 2 * accum[i];
    };

    let hashedKey = std::hash::pedersen(priv_key);
    [hashedKey]

}