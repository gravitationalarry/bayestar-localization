/*                                           !+:
                                            ~';j!
                                           '~~~t5~
                ';!!:~~''''.......        .~'~'jX;
            .>=+!!;;::~~~''''......       ~.'''J6~
          :cjJtc=>!;:~''''''''.......... ...'''Sc...
        'jtQQYYJc=>!;:~'''''.'''.........     >Y; ...
       ~tYNHMNQSJi+>>;:~''''.............''.. ci'...'.
      .jKMMNStcccci>:'''..............'~'.. ;i!.......
      !6MWKXSJci=>!;;~..'.........'~''... '=>'........
      iKKKWNY;''''';;'.''.                .=+'..........
    .cQMY~  .~'''+>''''  .    ......     .i>'...........
   =NX+.   '!=itKS'.''    .''''.... .    ~c!~..........'..
  !!''. .~!+6MNt~.''.    '''..          !i;:'..........'..
    '''~;JSc:''~'''. '.~.'...          =+!;~...............
    ':=6c!;:cMc.''. !':~.'::'.~'     ;c>!!;'...............
    tY!!!!!SWi'''. :; !'...'. .  . .j=>>!!:................
    'ici=+j5c'''. .:..~'.'.       +t>++>>!~.................
     'jYjiJj~'''  '~ .~'..'.    'Y=>++++>!'..''~'''.'..
       JYJ5!'''. .~' .'''..    =t!>>>+++>!'.:>!;:~~~'''.
      .+SSj''''  .~' .''.'   .Y+!!!>>++iYt~>>;;;::~~~'''''''..
       :QD>'''. ..'' ..'.   ;Y:;!!!!=6H5> ;+>!>>;'''''''''''''''''
      .'KX:'''  ..''  '    cc~::;!tNX+!!.;==j>.''''.           .''''.
      .'Q6~'''  ..''.'.   t!~::iKHi!>>+';=6!.''.  '~~::::~'.. ..   '''.
       .QJ'''' ...~.:.  'j~~>XWY>+==>+~'Jc '.  ~!+i+>!;;:~''..       .''
      ..Xt'''. .'' ~'  ~>;6MKc=iicii>;;6:'.  :cYSYtji>;;:~'''.   .     .'.
      . Y6'''. '  .'  ;5WW5ccjtJY56YiiY.'   =5WNNHXYc=>!:~'..'...       ..
       .iK~'''    .. =MQttJY5SDQHHH6JJ.'   >5WMMMNQ6j=!;:~'.......      ..
                 .' =tcJ5SDQHNNWWWNSY''  .:jXWMMNXti+;:~''.....'.      '
                 ~':>it6QHNWMMMMMWHD:.. :~;ttQWMWNHHXJi+!;:''........     ..
                 :'~;!=5KNWMMMMMWHQi.. ';'~cYQNWWNKD6ti+!;:~'..........   ..
                 . .'c5>jQNWWMWWNQ6:.  ;:.~+=jt5SS6Ytc=+!;:~'......'.   .
                 ..~;QHNM55KWWNNQ5='. .!:.'+=icjjjcci=+>!;:~'''''......
                  ''tNKHHXWH5+>=i=:.  .!:.'>>+iiii==+>!!;:~'..    .....
                  ''SHQHSXWHKKKQQ5i+;:~~~' ~+!>+>!!!;;:~. ....     ....
                  .'XHQKtQWNHKQQDXXD5t+;!~..~~:~~~~'....'''..'.   .''...
                  ''tWKX=HWNNHHQSYtc=!;~::'. ...     ...  .''''''''''''.
                  .  cWY>KNNHKD5j=>!;~~'.'..          ::'.  .~'.''''''''.
                  .   !X!5KKQSYc++!::''...             ::~.   .~'''''''''.

****************  ____  _____  ______________________    ____     **************
***************  / __ )/   \ \/ / ____/ ___/_  __/   |  / __ \   ***************
**************  / __  / /| |\  / __/  \__ \ / / / /| | / /_/ /  ****************
*************  / /_/ / ___ |/ / /___ ___/ // / / ___ |/ _, _/  *****************
************  /_____/_/  |_/_/_____//____//_/ /_/  |_/_/ |_|  ******************
*/


/*
 * Copyright (C) 2012  Leo Singer
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BAYESTAR_SKY_MAP_H
#define BAYESTAR_SKY_MAP_H


typedef enum
{
    BAYESTAR_PRIOR_UNIFORM_IN_LOG_DISTANCE,
    BAYESTAR_PRIOR_UNIFORM_IN_VOLUME
} bayestar_prior_t;


/* Perform sky localization based on TDOAs alone. */
double *bayestar_sky_map_tdoa(
    long *npix, /* In/out: number of HEALPix pixels. */
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    const double **locs, /* Input: array of detector positions. */
    const double *toas, /* Input: array of times of arrival. */
    const double *s2_toas /* Input: uncertainties in times of arrival. */
);

/* Perform sky localization based on TDOAs and amplitude. */
double *bayestar_sky_map_tdoa_snr(
    long *npix, /* In/out: number of HEALPix pixels. */
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    /* FIXME: make const; change XLALComputeDetAMResponse prototype */
    /* const */ float **responses, /* Pointers to detector responses. */
    const double **locations, /* Pointers to locations of detectors in Cartesian geographic coordinates. */
    const double *toas, /* Input: array of times of arrival with arbitrary relative offset. (Make toas[0] == 0.) */
    const double *snrs, /* Input: array of SNRs. */
    const double *s2_toas, /* Measurement variance of TOAs. */
    const double *horizons, /* Distances at which a source would produce an SNR of 1 in each detector. */
    double min_distance,
    double max_distance,
    bayestar_prior_t prior
);

#endif /* BAYESTAR_SKY_MAP_H */
