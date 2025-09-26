#include "ElectronRepulsionGeom0010ContrRecXXHI.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxhi(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxhi,
                                            const size_t idx_xxgi,
                                            const size_t idx_geom_10_xxgi,
                                            const size_t idx_geom_10_xxgk,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_cartesian_components(std::array<int, 1>{a_angmom,});

    const auto bcomps = tensor::number_of_cartesian_components(std::array<int, 1>{b_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : SSGI

            const auto gi_off = idx_xxgi + (i * bcomps + j) * 420;

            auto g_xxxx_xxxxxx = cbuffer.data(gi_off + 0);

            auto g_xxxx_xxxxxy = cbuffer.data(gi_off + 1);

            auto g_xxxx_xxxxxz = cbuffer.data(gi_off + 2);

            auto g_xxxx_xxxxyy = cbuffer.data(gi_off + 3);

            auto g_xxxx_xxxxyz = cbuffer.data(gi_off + 4);

            auto g_xxxx_xxxxzz = cbuffer.data(gi_off + 5);

            auto g_xxxx_xxxyyy = cbuffer.data(gi_off + 6);

            auto g_xxxx_xxxyyz = cbuffer.data(gi_off + 7);

            auto g_xxxx_xxxyzz = cbuffer.data(gi_off + 8);

            auto g_xxxx_xxxzzz = cbuffer.data(gi_off + 9);

            auto g_xxxx_xxyyyy = cbuffer.data(gi_off + 10);

            auto g_xxxx_xxyyyz = cbuffer.data(gi_off + 11);

            auto g_xxxx_xxyyzz = cbuffer.data(gi_off + 12);

            auto g_xxxx_xxyzzz = cbuffer.data(gi_off + 13);

            auto g_xxxx_xxzzzz = cbuffer.data(gi_off + 14);

            auto g_xxxx_xyyyyy = cbuffer.data(gi_off + 15);

            auto g_xxxx_xyyyyz = cbuffer.data(gi_off + 16);

            auto g_xxxx_xyyyzz = cbuffer.data(gi_off + 17);

            auto g_xxxx_xyyzzz = cbuffer.data(gi_off + 18);

            auto g_xxxx_xyzzzz = cbuffer.data(gi_off + 19);

            auto g_xxxx_xzzzzz = cbuffer.data(gi_off + 20);

            auto g_xxxx_yyyyyy = cbuffer.data(gi_off + 21);

            auto g_xxxx_yyyyyz = cbuffer.data(gi_off + 22);

            auto g_xxxx_yyyyzz = cbuffer.data(gi_off + 23);

            auto g_xxxx_yyyzzz = cbuffer.data(gi_off + 24);

            auto g_xxxx_yyzzzz = cbuffer.data(gi_off + 25);

            auto g_xxxx_yzzzzz = cbuffer.data(gi_off + 26);

            auto g_xxxx_zzzzzz = cbuffer.data(gi_off + 27);

            auto g_xyzz_xxxxyz = cbuffer.data(gi_off + 228);

            auto g_xyzz_xxxxzz = cbuffer.data(gi_off + 229);

            auto g_xyzz_xxxyyy = cbuffer.data(gi_off + 230);

            auto g_xyzz_xxxyyz = cbuffer.data(gi_off + 231);

            auto g_xyzz_xxxyzz = cbuffer.data(gi_off + 232);

            auto g_xyzz_xxxzzz = cbuffer.data(gi_off + 233);

            auto g_xyzz_xxyyyy = cbuffer.data(gi_off + 234);

            auto g_xyzz_xxyyyz = cbuffer.data(gi_off + 235);

            auto g_xyzz_xxyyzz = cbuffer.data(gi_off + 236);

            auto g_xyzz_xxyzzz = cbuffer.data(gi_off + 237);

            auto g_xyzz_xxzzzz = cbuffer.data(gi_off + 238);

            auto g_xyzz_xyyyyy = cbuffer.data(gi_off + 239);

            auto g_xyzz_xyyyyz = cbuffer.data(gi_off + 240);

            auto g_xyzz_xyyyzz = cbuffer.data(gi_off + 241);

            auto g_xyzz_xyyzzz = cbuffer.data(gi_off + 242);

            auto g_xyzz_xyzzzz = cbuffer.data(gi_off + 243);

            auto g_xyzz_xzzzzz = cbuffer.data(gi_off + 244);

            auto g_xyzz_yyyyyy = cbuffer.data(gi_off + 245);

            auto g_xyzz_yyyyyz = cbuffer.data(gi_off + 246);

            auto g_xyzz_yyyyzz = cbuffer.data(gi_off + 247);

            auto g_xyzz_yyyzzz = cbuffer.data(gi_off + 248);

            auto g_xyzz_yyzzzz = cbuffer.data(gi_off + 249);

            auto g_xyzz_yzzzzz = cbuffer.data(gi_off + 250);

            auto g_xyzz_zzzzzz = cbuffer.data(gi_off + 251);

            auto g_xzzz_xxxxxx = cbuffer.data(gi_off + 252);

            auto g_xzzz_xxxxxy = cbuffer.data(gi_off + 253);

            auto g_xzzz_xxxxxz = cbuffer.data(gi_off + 254);

            auto g_xzzz_xxxxyy = cbuffer.data(gi_off + 255);

            auto g_xzzz_xxxxyz = cbuffer.data(gi_off + 256);

            auto g_xzzz_xxxxzz = cbuffer.data(gi_off + 257);

            auto g_xzzz_xxxyyy = cbuffer.data(gi_off + 258);

            auto g_xzzz_xxxyyz = cbuffer.data(gi_off + 259);

            auto g_xzzz_xxxyzz = cbuffer.data(gi_off + 260);

            auto g_xzzz_xxxzzz = cbuffer.data(gi_off + 261);

            auto g_xzzz_xxyyyy = cbuffer.data(gi_off + 262);

            auto g_xzzz_xxyyyz = cbuffer.data(gi_off + 263);

            auto g_xzzz_xxyyzz = cbuffer.data(gi_off + 264);

            auto g_xzzz_xxyzzz = cbuffer.data(gi_off + 265);

            auto g_xzzz_xxzzzz = cbuffer.data(gi_off + 266);

            auto g_xzzz_xyyyyy = cbuffer.data(gi_off + 267);

            auto g_xzzz_xyyyyz = cbuffer.data(gi_off + 268);

            auto g_xzzz_xyyyzz = cbuffer.data(gi_off + 269);

            auto g_xzzz_xyyzzz = cbuffer.data(gi_off + 270);

            auto g_xzzz_xyzzzz = cbuffer.data(gi_off + 271);

            auto g_xzzz_xzzzzz = cbuffer.data(gi_off + 272);

            auto g_xzzz_yyyyyy = cbuffer.data(gi_off + 273);

            auto g_xzzz_yyyyyz = cbuffer.data(gi_off + 274);

            auto g_xzzz_yyyyzz = cbuffer.data(gi_off + 275);

            auto g_xzzz_yyyzzz = cbuffer.data(gi_off + 276);

            auto g_xzzz_yyzzzz = cbuffer.data(gi_off + 277);

            auto g_xzzz_yzzzzz = cbuffer.data(gi_off + 278);

            auto g_xzzz_zzzzzz = cbuffer.data(gi_off + 279);

            auto g_yyyy_xxxxxx = cbuffer.data(gi_off + 280);

            auto g_yyyy_xxxxxy = cbuffer.data(gi_off + 281);

            auto g_yyyy_xxxxxz = cbuffer.data(gi_off + 282);

            auto g_yyyy_xxxxyy = cbuffer.data(gi_off + 283);

            auto g_yyyy_xxxxyz = cbuffer.data(gi_off + 284);

            auto g_yyyy_xxxxzz = cbuffer.data(gi_off + 285);

            auto g_yyyy_xxxyyy = cbuffer.data(gi_off + 286);

            auto g_yyyy_xxxyyz = cbuffer.data(gi_off + 287);

            auto g_yyyy_xxxyzz = cbuffer.data(gi_off + 288);

            auto g_yyyy_xxxzzz = cbuffer.data(gi_off + 289);

            auto g_yyyy_xxyyyy = cbuffer.data(gi_off + 290);

            auto g_yyyy_xxyyyz = cbuffer.data(gi_off + 291);

            auto g_yyyy_xxyyzz = cbuffer.data(gi_off + 292);

            auto g_yyyy_xxyzzz = cbuffer.data(gi_off + 293);

            auto g_yyyy_xxzzzz = cbuffer.data(gi_off + 294);

            auto g_yyyy_xyyyyy = cbuffer.data(gi_off + 295);

            auto g_yyyy_xyyyyz = cbuffer.data(gi_off + 296);

            auto g_yyyy_xyyyzz = cbuffer.data(gi_off + 297);

            auto g_yyyy_xyyzzz = cbuffer.data(gi_off + 298);

            auto g_yyyy_xyzzzz = cbuffer.data(gi_off + 299);

            auto g_yyyy_xzzzzz = cbuffer.data(gi_off + 300);

            auto g_yyyy_yyyyyy = cbuffer.data(gi_off + 301);

            auto g_yyyy_yyyyyz = cbuffer.data(gi_off + 302);

            auto g_yyyy_yyyyzz = cbuffer.data(gi_off + 303);

            auto g_yyyy_yyyzzz = cbuffer.data(gi_off + 304);

            auto g_yyyy_yyzzzz = cbuffer.data(gi_off + 305);

            auto g_yyyy_yzzzzz = cbuffer.data(gi_off + 306);

            auto g_yyyy_zzzzzz = cbuffer.data(gi_off + 307);

            auto g_yyyz_xxxxxx = cbuffer.data(gi_off + 308);

            auto g_yyyz_xxxxxy = cbuffer.data(gi_off + 309);

            auto g_yyyz_xxxxxz = cbuffer.data(gi_off + 310);

            auto g_yyyz_xxxxyy = cbuffer.data(gi_off + 311);

            auto g_yyyz_xxxxyz = cbuffer.data(gi_off + 312);

            auto g_yyyz_xxxxzz = cbuffer.data(gi_off + 313);

            auto g_yyyz_xxxyyy = cbuffer.data(gi_off + 314);

            auto g_yyyz_xxxyyz = cbuffer.data(gi_off + 315);

            auto g_yyyz_xxxyzz = cbuffer.data(gi_off + 316);

            auto g_yyyz_xxxzzz = cbuffer.data(gi_off + 317);

            auto g_yyyz_xxyyyy = cbuffer.data(gi_off + 318);

            auto g_yyyz_xxyyyz = cbuffer.data(gi_off + 319);

            auto g_yyyz_xxyyzz = cbuffer.data(gi_off + 320);

            auto g_yyyz_xxyzzz = cbuffer.data(gi_off + 321);

            auto g_yyyz_xxzzzz = cbuffer.data(gi_off + 322);

            auto g_yyyz_xyyyyy = cbuffer.data(gi_off + 323);

            auto g_yyyz_xyyyyz = cbuffer.data(gi_off + 324);

            auto g_yyyz_xyyyzz = cbuffer.data(gi_off + 325);

            auto g_yyyz_xyyzzz = cbuffer.data(gi_off + 326);

            auto g_yyyz_xyzzzz = cbuffer.data(gi_off + 327);

            auto g_yyyz_xzzzzz = cbuffer.data(gi_off + 328);

            auto g_yyyz_yyyyyy = cbuffer.data(gi_off + 329);

            auto g_yyyz_yyyyyz = cbuffer.data(gi_off + 330);

            auto g_yyyz_yyyyzz = cbuffer.data(gi_off + 331);

            auto g_yyyz_yyyzzz = cbuffer.data(gi_off + 332);

            auto g_yyyz_yyzzzz = cbuffer.data(gi_off + 333);

            auto g_yyyz_yzzzzz = cbuffer.data(gi_off + 334);

            auto g_yyyz_zzzzzz = cbuffer.data(gi_off + 335);

            auto g_yyzz_xxxxxx = cbuffer.data(gi_off + 336);

            auto g_yyzz_xxxxxy = cbuffer.data(gi_off + 337);

            auto g_yyzz_xxxxxz = cbuffer.data(gi_off + 338);

            auto g_yyzz_xxxxyy = cbuffer.data(gi_off + 339);

            auto g_yyzz_xxxxyz = cbuffer.data(gi_off + 340);

            auto g_yyzz_xxxxzz = cbuffer.data(gi_off + 341);

            auto g_yyzz_xxxyyy = cbuffer.data(gi_off + 342);

            auto g_yyzz_xxxyyz = cbuffer.data(gi_off + 343);

            auto g_yyzz_xxxyzz = cbuffer.data(gi_off + 344);

            auto g_yyzz_xxxzzz = cbuffer.data(gi_off + 345);

            auto g_yyzz_xxyyyy = cbuffer.data(gi_off + 346);

            auto g_yyzz_xxyyyz = cbuffer.data(gi_off + 347);

            auto g_yyzz_xxyyzz = cbuffer.data(gi_off + 348);

            auto g_yyzz_xxyzzz = cbuffer.data(gi_off + 349);

            auto g_yyzz_xxzzzz = cbuffer.data(gi_off + 350);

            auto g_yyzz_xyyyyy = cbuffer.data(gi_off + 351);

            auto g_yyzz_xyyyyz = cbuffer.data(gi_off + 352);

            auto g_yyzz_xyyyzz = cbuffer.data(gi_off + 353);

            auto g_yyzz_xyyzzz = cbuffer.data(gi_off + 354);

            auto g_yyzz_xyzzzz = cbuffer.data(gi_off + 355);

            auto g_yyzz_xzzzzz = cbuffer.data(gi_off + 356);

            auto g_yyzz_yyyyyy = cbuffer.data(gi_off + 357);

            auto g_yyzz_yyyyyz = cbuffer.data(gi_off + 358);

            auto g_yyzz_yyyyzz = cbuffer.data(gi_off + 359);

            auto g_yyzz_yyyzzz = cbuffer.data(gi_off + 360);

            auto g_yyzz_yyzzzz = cbuffer.data(gi_off + 361);

            auto g_yyzz_yzzzzz = cbuffer.data(gi_off + 362);

            auto g_yyzz_zzzzzz = cbuffer.data(gi_off + 363);

            auto g_yzzz_xxxxxx = cbuffer.data(gi_off + 364);

            auto g_yzzz_xxxxxy = cbuffer.data(gi_off + 365);

            auto g_yzzz_xxxxxz = cbuffer.data(gi_off + 366);

            auto g_yzzz_xxxxyy = cbuffer.data(gi_off + 367);

            auto g_yzzz_xxxxyz = cbuffer.data(gi_off + 368);

            auto g_yzzz_xxxxzz = cbuffer.data(gi_off + 369);

            auto g_yzzz_xxxyyy = cbuffer.data(gi_off + 370);

            auto g_yzzz_xxxyyz = cbuffer.data(gi_off + 371);

            auto g_yzzz_xxxyzz = cbuffer.data(gi_off + 372);

            auto g_yzzz_xxxzzz = cbuffer.data(gi_off + 373);

            auto g_yzzz_xxyyyy = cbuffer.data(gi_off + 374);

            auto g_yzzz_xxyyyz = cbuffer.data(gi_off + 375);

            auto g_yzzz_xxyyzz = cbuffer.data(gi_off + 376);

            auto g_yzzz_xxyzzz = cbuffer.data(gi_off + 377);

            auto g_yzzz_xxzzzz = cbuffer.data(gi_off + 378);

            auto g_yzzz_xyyyyy = cbuffer.data(gi_off + 379);

            auto g_yzzz_xyyyyz = cbuffer.data(gi_off + 380);

            auto g_yzzz_xyyyzz = cbuffer.data(gi_off + 381);

            auto g_yzzz_xyyzzz = cbuffer.data(gi_off + 382);

            auto g_yzzz_xyzzzz = cbuffer.data(gi_off + 383);

            auto g_yzzz_xzzzzz = cbuffer.data(gi_off + 384);

            auto g_yzzz_yyyyyy = cbuffer.data(gi_off + 385);

            auto g_yzzz_yyyyyz = cbuffer.data(gi_off + 386);

            auto g_yzzz_yyyyzz = cbuffer.data(gi_off + 387);

            auto g_yzzz_yyyzzz = cbuffer.data(gi_off + 388);

            auto g_yzzz_yyzzzz = cbuffer.data(gi_off + 389);

            auto g_yzzz_yzzzzz = cbuffer.data(gi_off + 390);

            auto g_yzzz_zzzzzz = cbuffer.data(gi_off + 391);

            auto g_zzzz_xxxxxx = cbuffer.data(gi_off + 392);

            auto g_zzzz_xxxxxy = cbuffer.data(gi_off + 393);

            auto g_zzzz_xxxxxz = cbuffer.data(gi_off + 394);

            auto g_zzzz_xxxxyy = cbuffer.data(gi_off + 395);

            auto g_zzzz_xxxxyz = cbuffer.data(gi_off + 396);

            auto g_zzzz_xxxxzz = cbuffer.data(gi_off + 397);

            auto g_zzzz_xxxyyy = cbuffer.data(gi_off + 398);

            auto g_zzzz_xxxyyz = cbuffer.data(gi_off + 399);

            auto g_zzzz_xxxyzz = cbuffer.data(gi_off + 400);

            auto g_zzzz_xxxzzz = cbuffer.data(gi_off + 401);

            auto g_zzzz_xxyyyy = cbuffer.data(gi_off + 402);

            auto g_zzzz_xxyyyz = cbuffer.data(gi_off + 403);

            auto g_zzzz_xxyyzz = cbuffer.data(gi_off + 404);

            auto g_zzzz_xxyzzz = cbuffer.data(gi_off + 405);

            auto g_zzzz_xxzzzz = cbuffer.data(gi_off + 406);

            auto g_zzzz_xyyyyy = cbuffer.data(gi_off + 407);

            auto g_zzzz_xyyyyz = cbuffer.data(gi_off + 408);

            auto g_zzzz_xyyyzz = cbuffer.data(gi_off + 409);

            auto g_zzzz_xyyzzz = cbuffer.data(gi_off + 410);

            auto g_zzzz_xyzzzz = cbuffer.data(gi_off + 411);

            auto g_zzzz_xzzzzz = cbuffer.data(gi_off + 412);

            auto g_zzzz_yyyyyy = cbuffer.data(gi_off + 413);

            auto g_zzzz_yyyyyz = cbuffer.data(gi_off + 414);

            auto g_zzzz_yyyyzz = cbuffer.data(gi_off + 415);

            auto g_zzzz_yyyzzz = cbuffer.data(gi_off + 416);

            auto g_zzzz_yyzzzz = cbuffer.data(gi_off + 417);

            auto g_zzzz_yzzzzz = cbuffer.data(gi_off + 418);

            auto g_zzzz_zzzzzz = cbuffer.data(gi_off + 419);

            /// Set up components of auxilary buffer : SSGI

            const auto gi_geom_10_off = idx_geom_10_xxgi + (i * bcomps + j) * 420;

            auto g_x_0_xxxx_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxx_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxx_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxx_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxx_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxx_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxx_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxx_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxx_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxx_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxx_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxx_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxx_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxx_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxx_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxx_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxx_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxx_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxx_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxx_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxx_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxx_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxx_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxx_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxx_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxx_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxx_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxx_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxy_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxy_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxxy_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxy_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxy_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxy_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxy_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxy_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxy_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxy_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxy_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxy_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxy_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxy_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxxy_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxy_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxy_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxxy_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxy_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxy_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxxy_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxy_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxxy_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxy_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxy_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxy_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxy_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxy_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxxz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xxxz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxxz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxxz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxxz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxxz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxxz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxxz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxxz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxxz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxxz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxxz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxxz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxxz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxxz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxxz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxxz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxxz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxxz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxxz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxxz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxxz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxxz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxxz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxxz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xxyy_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxyy_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxyy_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxyy_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxyy_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxyy_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xxyy_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxyy_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxyy_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxyy_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxyy_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxyy_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxyy_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxyy_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxyy_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxyy_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xxyy_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xxyy_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xxyy_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xxyy_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xxyy_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xxyy_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xxyy_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xxyy_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xxyy_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xxyy_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xxyy_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xxyy_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xxyz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xxyz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xxyz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xxyz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xxyz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xxyz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xxyz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xxyz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xxyz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xxyz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xxyz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xxyz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xxyz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xxyz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xxyz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xxyz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xxyz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xxyz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xxyz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xxyz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xxyz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xxyz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xxyz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xxyz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xxyz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xxyz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xxyz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xxyz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_xxzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xxzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xxzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xxzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xxzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xxzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xxzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_xxzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xxzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xxzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_xxzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xxzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xxzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_xxzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xxzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xxzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xxzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_xxzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xxzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xxzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_xxzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xxzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_xxzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xxzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xxzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_xxzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xxzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xxzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_x_0_xyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_xyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_xyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_xyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_xyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_xyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_xyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_xyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_xyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_xyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_xyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_xyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_xyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_xyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_xyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_xyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_xyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_xyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_xyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_xyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_xyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_xyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_xyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_xyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_xyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_xyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_xyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_xyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_xyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_xyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_xyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_xyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_xyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_xyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_xyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_xyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_xyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_xyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_xyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_xyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_xyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_xyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 209);

            auto g_x_0_xyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_xyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_xyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_xyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_xyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_xyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_xyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_xyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_xyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_xyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_xyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_xyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_xyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_xyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 223);

            auto g_x_0_xyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 224);

            auto g_x_0_xyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_xyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_xyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_xyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_xyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 229);

            auto g_x_0_xyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 230);

            auto g_x_0_xyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_xyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_xyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 233);

            auto g_x_0_xyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_xyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_xyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_xyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_xyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_xyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 239);

            auto g_x_0_xyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_xyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_xyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_xyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_xyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_xyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_xyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_xyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_xyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_xyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_xyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_xyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 251);

            auto g_x_0_xzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 252);

            auto g_x_0_xzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_xzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 254);

            auto g_x_0_xzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_xzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_xzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 257);

            auto g_x_0_xzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_xzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 259);

            auto g_x_0_xzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_xzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 261);

            auto g_x_0_xzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_xzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_xzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_xzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_xzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 266);

            auto g_x_0_xzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_xzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_xzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 269);

            auto g_x_0_xzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_xzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_xzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 272);

            auto g_x_0_xzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 273);

            auto g_x_0_xzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_xzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 275);

            auto g_x_0_xzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_xzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_xzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 278);

            auto g_x_0_xzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 279);

            auto g_x_0_yyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 280);

            auto g_x_0_yyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 281);

            auto g_x_0_yyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 282);

            auto g_x_0_yyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 283);

            auto g_x_0_yyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 284);

            auto g_x_0_yyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 285);

            auto g_x_0_yyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 286);

            auto g_x_0_yyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 287);

            auto g_x_0_yyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 288);

            auto g_x_0_yyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 289);

            auto g_x_0_yyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 290);

            auto g_x_0_yyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 291);

            auto g_x_0_yyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 292);

            auto g_x_0_yyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 293);

            auto g_x_0_yyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 294);

            auto g_x_0_yyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 295);

            auto g_x_0_yyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 296);

            auto g_x_0_yyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 297);

            auto g_x_0_yyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 298);

            auto g_x_0_yyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 299);

            auto g_x_0_yyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 300);

            auto g_x_0_yyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 301);

            auto g_x_0_yyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 302);

            auto g_x_0_yyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 303);

            auto g_x_0_yyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 304);

            auto g_x_0_yyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 305);

            auto g_x_0_yyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 306);

            auto g_x_0_yyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 307);

            auto g_x_0_yyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 308);

            auto g_x_0_yyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 309);

            auto g_x_0_yyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 310);

            auto g_x_0_yyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 311);

            auto g_x_0_yyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 312);

            auto g_x_0_yyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 313);

            auto g_x_0_yyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 314);

            auto g_x_0_yyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 315);

            auto g_x_0_yyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 316);

            auto g_x_0_yyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 317);

            auto g_x_0_yyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 318);

            auto g_x_0_yyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 319);

            auto g_x_0_yyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 320);

            auto g_x_0_yyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 321);

            auto g_x_0_yyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 322);

            auto g_x_0_yyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 323);

            auto g_x_0_yyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 324);

            auto g_x_0_yyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 325);

            auto g_x_0_yyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 326);

            auto g_x_0_yyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 327);

            auto g_x_0_yyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 328);

            auto g_x_0_yyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 329);

            auto g_x_0_yyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 330);

            auto g_x_0_yyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 331);

            auto g_x_0_yyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 332);

            auto g_x_0_yyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 333);

            auto g_x_0_yyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 334);

            auto g_x_0_yyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 335);

            auto g_x_0_yyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 336);

            auto g_x_0_yyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 337);

            auto g_x_0_yyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 338);

            auto g_x_0_yyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 339);

            auto g_x_0_yyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 340);

            auto g_x_0_yyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 341);

            auto g_x_0_yyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 342);

            auto g_x_0_yyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 343);

            auto g_x_0_yyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 344);

            auto g_x_0_yyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 345);

            auto g_x_0_yyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 346);

            auto g_x_0_yyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 347);

            auto g_x_0_yyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 348);

            auto g_x_0_yyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 349);

            auto g_x_0_yyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 350);

            auto g_x_0_yyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 351);

            auto g_x_0_yyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 352);

            auto g_x_0_yyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 353);

            auto g_x_0_yyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 354);

            auto g_x_0_yyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 355);

            auto g_x_0_yyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 356);

            auto g_x_0_yyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 357);

            auto g_x_0_yyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 358);

            auto g_x_0_yyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 359);

            auto g_x_0_yyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 360);

            auto g_x_0_yyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 361);

            auto g_x_0_yyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 362);

            auto g_x_0_yyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 363);

            auto g_x_0_yzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 364);

            auto g_x_0_yzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 365);

            auto g_x_0_yzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 366);

            auto g_x_0_yzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 367);

            auto g_x_0_yzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 368);

            auto g_x_0_yzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 369);

            auto g_x_0_yzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 370);

            auto g_x_0_yzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 371);

            auto g_x_0_yzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 372);

            auto g_x_0_yzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 373);

            auto g_x_0_yzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 374);

            auto g_x_0_yzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 375);

            auto g_x_0_yzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 376);

            auto g_x_0_yzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 377);

            auto g_x_0_yzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 378);

            auto g_x_0_yzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 379);

            auto g_x_0_yzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 380);

            auto g_x_0_yzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 381);

            auto g_x_0_yzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 382);

            auto g_x_0_yzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 383);

            auto g_x_0_yzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 384);

            auto g_x_0_yzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 385);

            auto g_x_0_yzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 386);

            auto g_x_0_yzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 387);

            auto g_x_0_yzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 388);

            auto g_x_0_yzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 389);

            auto g_x_0_yzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 390);

            auto g_x_0_yzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 391);

            auto g_x_0_zzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 392);

            auto g_x_0_zzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 393);

            auto g_x_0_zzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 394);

            auto g_x_0_zzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 395);

            auto g_x_0_zzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 396);

            auto g_x_0_zzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 397);

            auto g_x_0_zzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 398);

            auto g_x_0_zzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 399);

            auto g_x_0_zzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 400);

            auto g_x_0_zzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 401);

            auto g_x_0_zzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 402);

            auto g_x_0_zzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 403);

            auto g_x_0_zzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 404);

            auto g_x_0_zzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 405);

            auto g_x_0_zzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 406);

            auto g_x_0_zzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 407);

            auto g_x_0_zzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 408);

            auto g_x_0_zzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 409);

            auto g_x_0_zzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 410);

            auto g_x_0_zzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 411);

            auto g_x_0_zzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 412);

            auto g_x_0_zzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 413);

            auto g_x_0_zzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 414);

            auto g_x_0_zzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 415);

            auto g_x_0_zzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 416);

            auto g_x_0_zzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 417);

            auto g_x_0_zzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 418);

            auto g_x_0_zzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 419);

            auto g_y_0_xxxx_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 0);

            auto g_y_0_xxxx_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 1);

            auto g_y_0_xxxx_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 2);

            auto g_y_0_xxxx_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 3);

            auto g_y_0_xxxx_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 4);

            auto g_y_0_xxxx_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 5);

            auto g_y_0_xxxx_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 6);

            auto g_y_0_xxxx_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 7);

            auto g_y_0_xxxx_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 8);

            auto g_y_0_xxxx_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 9);

            auto g_y_0_xxxx_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 10);

            auto g_y_0_xxxx_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 11);

            auto g_y_0_xxxx_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 12);

            auto g_y_0_xxxx_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 13);

            auto g_y_0_xxxx_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 14);

            auto g_y_0_xxxx_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 15);

            auto g_y_0_xxxx_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 16);

            auto g_y_0_xxxx_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 17);

            auto g_y_0_xxxx_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 18);

            auto g_y_0_xxxx_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 19);

            auto g_y_0_xxxx_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 20);

            auto g_y_0_xxxx_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 21);

            auto g_y_0_xxxx_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 22);

            auto g_y_0_xxxx_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 23);

            auto g_y_0_xxxx_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 24);

            auto g_y_0_xxxx_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 25);

            auto g_y_0_xxxx_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 26);

            auto g_y_0_xxxx_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 27);

            auto g_y_0_xxxy_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 28);

            auto g_y_0_xxxy_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 29);

            auto g_y_0_xxxy_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 30);

            auto g_y_0_xxxy_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 31);

            auto g_y_0_xxxy_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 32);

            auto g_y_0_xxxy_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 33);

            auto g_y_0_xxxy_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 34);

            auto g_y_0_xxxy_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 35);

            auto g_y_0_xxxy_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 36);

            auto g_y_0_xxxy_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 37);

            auto g_y_0_xxxy_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 38);

            auto g_y_0_xxxy_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 39);

            auto g_y_0_xxxy_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 40);

            auto g_y_0_xxxy_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 41);

            auto g_y_0_xxxy_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 42);

            auto g_y_0_xxxy_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 43);

            auto g_y_0_xxxy_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 44);

            auto g_y_0_xxxy_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 45);

            auto g_y_0_xxxy_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 46);

            auto g_y_0_xxxy_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 47);

            auto g_y_0_xxxy_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 48);

            auto g_y_0_xxxy_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 49);

            auto g_y_0_xxxy_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 50);

            auto g_y_0_xxxy_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 51);

            auto g_y_0_xxxy_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 52);

            auto g_y_0_xxxy_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 53);

            auto g_y_0_xxxy_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 54);

            auto g_y_0_xxxy_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 55);

            auto g_y_0_xxxz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 56);

            auto g_y_0_xxxz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 57);

            auto g_y_0_xxxz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 58);

            auto g_y_0_xxxz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 59);

            auto g_y_0_xxxz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 60);

            auto g_y_0_xxxz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 61);

            auto g_y_0_xxxz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 62);

            auto g_y_0_xxxz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 63);

            auto g_y_0_xxxz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 64);

            auto g_y_0_xxxz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 65);

            auto g_y_0_xxxz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 66);

            auto g_y_0_xxxz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 67);

            auto g_y_0_xxxz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 68);

            auto g_y_0_xxxz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 69);

            auto g_y_0_xxxz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 70);

            auto g_y_0_xxxz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 71);

            auto g_y_0_xxxz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 72);

            auto g_y_0_xxxz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 73);

            auto g_y_0_xxxz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 74);

            auto g_y_0_xxxz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 75);

            auto g_y_0_xxxz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 76);

            auto g_y_0_xxxz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 77);

            auto g_y_0_xxxz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 78);

            auto g_y_0_xxxz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 79);

            auto g_y_0_xxxz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 80);

            auto g_y_0_xxxz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 81);

            auto g_y_0_xxxz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 82);

            auto g_y_0_xxxz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 83);

            auto g_y_0_xxyy_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 84);

            auto g_y_0_xxyy_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 85);

            auto g_y_0_xxyy_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 86);

            auto g_y_0_xxyy_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 87);

            auto g_y_0_xxyy_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 88);

            auto g_y_0_xxyy_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 89);

            auto g_y_0_xxyy_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 90);

            auto g_y_0_xxyy_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 91);

            auto g_y_0_xxyy_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 92);

            auto g_y_0_xxyy_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 93);

            auto g_y_0_xxyy_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 94);

            auto g_y_0_xxyy_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 95);

            auto g_y_0_xxyy_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 96);

            auto g_y_0_xxyy_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 97);

            auto g_y_0_xxyy_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 98);

            auto g_y_0_xxyy_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 99);

            auto g_y_0_xxyy_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 100);

            auto g_y_0_xxyy_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 101);

            auto g_y_0_xxyy_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 102);

            auto g_y_0_xxyy_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 103);

            auto g_y_0_xxyy_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 104);

            auto g_y_0_xxyy_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 105);

            auto g_y_0_xxyy_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 106);

            auto g_y_0_xxyy_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 107);

            auto g_y_0_xxyy_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 108);

            auto g_y_0_xxyy_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 109);

            auto g_y_0_xxyy_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 110);

            auto g_y_0_xxyy_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 111);

            auto g_y_0_xxyz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 112);

            auto g_y_0_xxyz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 113);

            auto g_y_0_xxyz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 114);

            auto g_y_0_xxyz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 115);

            auto g_y_0_xxyz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 116);

            auto g_y_0_xxyz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 117);

            auto g_y_0_xxyz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 118);

            auto g_y_0_xxyz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 119);

            auto g_y_0_xxyz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 120);

            auto g_y_0_xxyz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 121);

            auto g_y_0_xxyz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 122);

            auto g_y_0_xxyz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 123);

            auto g_y_0_xxyz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 124);

            auto g_y_0_xxyz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 125);

            auto g_y_0_xxyz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 126);

            auto g_y_0_xxyz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 127);

            auto g_y_0_xxyz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 128);

            auto g_y_0_xxyz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 129);

            auto g_y_0_xxyz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 130);

            auto g_y_0_xxyz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 131);

            auto g_y_0_xxyz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 132);

            auto g_y_0_xxyz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 133);

            auto g_y_0_xxyz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 134);

            auto g_y_0_xxyz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 135);

            auto g_y_0_xxyz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 136);

            auto g_y_0_xxyz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 137);

            auto g_y_0_xxyz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 138);

            auto g_y_0_xxyz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 139);

            auto g_y_0_xxzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 140);

            auto g_y_0_xxzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 141);

            auto g_y_0_xxzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 142);

            auto g_y_0_xxzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 143);

            auto g_y_0_xxzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 144);

            auto g_y_0_xxzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 145);

            auto g_y_0_xxzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 146);

            auto g_y_0_xxzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 147);

            auto g_y_0_xxzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 148);

            auto g_y_0_xxzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 149);

            auto g_y_0_xxzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 150);

            auto g_y_0_xxzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 151);

            auto g_y_0_xxzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 152);

            auto g_y_0_xxzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 153);

            auto g_y_0_xxzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 154);

            auto g_y_0_xxzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 155);

            auto g_y_0_xxzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 156);

            auto g_y_0_xxzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 157);

            auto g_y_0_xxzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 158);

            auto g_y_0_xxzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 159);

            auto g_y_0_xxzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 160);

            auto g_y_0_xxzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 161);

            auto g_y_0_xxzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 162);

            auto g_y_0_xxzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 163);

            auto g_y_0_xxzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 164);

            auto g_y_0_xxzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 165);

            auto g_y_0_xxzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 166);

            auto g_y_0_xxzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 167);

            auto g_y_0_xyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 168);

            auto g_y_0_xyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 169);

            auto g_y_0_xyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 170);

            auto g_y_0_xyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 171);

            auto g_y_0_xyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 172);

            auto g_y_0_xyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 173);

            auto g_y_0_xyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 174);

            auto g_y_0_xyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 175);

            auto g_y_0_xyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 176);

            auto g_y_0_xyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 177);

            auto g_y_0_xyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 178);

            auto g_y_0_xyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 179);

            auto g_y_0_xyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 180);

            auto g_y_0_xyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 181);

            auto g_y_0_xyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 182);

            auto g_y_0_xyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 183);

            auto g_y_0_xyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 184);

            auto g_y_0_xyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 185);

            auto g_y_0_xyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 186);

            auto g_y_0_xyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 187);

            auto g_y_0_xyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 188);

            auto g_y_0_xyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 189);

            auto g_y_0_xyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 190);

            auto g_y_0_xyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 191);

            auto g_y_0_xyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 192);

            auto g_y_0_xyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 193);

            auto g_y_0_xyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 194);

            auto g_y_0_xyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 195);

            auto g_y_0_xyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 196);

            auto g_y_0_xyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 197);

            auto g_y_0_xyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 198);

            auto g_y_0_xyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 199);

            auto g_y_0_xyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 200);

            auto g_y_0_xyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 201);

            auto g_y_0_xyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 202);

            auto g_y_0_xyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 203);

            auto g_y_0_xyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 204);

            auto g_y_0_xyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 205);

            auto g_y_0_xyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 206);

            auto g_y_0_xyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 207);

            auto g_y_0_xyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 208);

            auto g_y_0_xyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 209);

            auto g_y_0_xyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 210);

            auto g_y_0_xyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 211);

            auto g_y_0_xyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 212);

            auto g_y_0_xyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 213);

            auto g_y_0_xyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 214);

            auto g_y_0_xyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 215);

            auto g_y_0_xyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 216);

            auto g_y_0_xyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 217);

            auto g_y_0_xyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 218);

            auto g_y_0_xyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 219);

            auto g_y_0_xyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 220);

            auto g_y_0_xyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 221);

            auto g_y_0_xyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 222);

            auto g_y_0_xyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 223);

            auto g_y_0_xyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 224);

            auto g_y_0_xyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 225);

            auto g_y_0_xyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 226);

            auto g_y_0_xyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 227);

            auto g_y_0_xyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 228);

            auto g_y_0_xyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 229);

            auto g_y_0_xyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 230);

            auto g_y_0_xyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 231);

            auto g_y_0_xyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 232);

            auto g_y_0_xyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 233);

            auto g_y_0_xyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 234);

            auto g_y_0_xyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 235);

            auto g_y_0_xyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 236);

            auto g_y_0_xyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 237);

            auto g_y_0_xyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 238);

            auto g_y_0_xyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 239);

            auto g_y_0_xyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 240);

            auto g_y_0_xyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 241);

            auto g_y_0_xyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 242);

            auto g_y_0_xyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 243);

            auto g_y_0_xyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 244);

            auto g_y_0_xyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 245);

            auto g_y_0_xyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 246);

            auto g_y_0_xyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 247);

            auto g_y_0_xyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 248);

            auto g_y_0_xyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 249);

            auto g_y_0_xyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 250);

            auto g_y_0_xyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 251);

            auto g_y_0_xzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 252);

            auto g_y_0_xzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 253);

            auto g_y_0_xzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 254);

            auto g_y_0_xzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 255);

            auto g_y_0_xzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 256);

            auto g_y_0_xzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 257);

            auto g_y_0_xzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 258);

            auto g_y_0_xzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 259);

            auto g_y_0_xzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 260);

            auto g_y_0_xzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 261);

            auto g_y_0_xzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 262);

            auto g_y_0_xzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 263);

            auto g_y_0_xzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 264);

            auto g_y_0_xzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 265);

            auto g_y_0_xzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 266);

            auto g_y_0_xzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 267);

            auto g_y_0_xzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 268);

            auto g_y_0_xzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 269);

            auto g_y_0_xzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 270);

            auto g_y_0_xzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 271);

            auto g_y_0_xzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 272);

            auto g_y_0_xzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 273);

            auto g_y_0_xzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 274);

            auto g_y_0_xzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 275);

            auto g_y_0_xzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 276);

            auto g_y_0_xzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 277);

            auto g_y_0_xzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 278);

            auto g_y_0_xzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 279);

            auto g_y_0_yyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 280);

            auto g_y_0_yyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 281);

            auto g_y_0_yyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 282);

            auto g_y_0_yyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 283);

            auto g_y_0_yyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 284);

            auto g_y_0_yyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 285);

            auto g_y_0_yyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 286);

            auto g_y_0_yyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 287);

            auto g_y_0_yyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 288);

            auto g_y_0_yyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 289);

            auto g_y_0_yyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 290);

            auto g_y_0_yyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 291);

            auto g_y_0_yyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 292);

            auto g_y_0_yyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 293);

            auto g_y_0_yyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 294);

            auto g_y_0_yyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 295);

            auto g_y_0_yyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 296);

            auto g_y_0_yyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 297);

            auto g_y_0_yyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 298);

            auto g_y_0_yyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 299);

            auto g_y_0_yyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 300);

            auto g_y_0_yyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 301);

            auto g_y_0_yyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 302);

            auto g_y_0_yyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 303);

            auto g_y_0_yyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 304);

            auto g_y_0_yyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 305);

            auto g_y_0_yyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 306);

            auto g_y_0_yyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 307);

            auto g_y_0_yyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 308);

            auto g_y_0_yyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 309);

            auto g_y_0_yyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 310);

            auto g_y_0_yyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 311);

            auto g_y_0_yyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 312);

            auto g_y_0_yyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 313);

            auto g_y_0_yyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 314);

            auto g_y_0_yyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 315);

            auto g_y_0_yyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 316);

            auto g_y_0_yyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 317);

            auto g_y_0_yyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 318);

            auto g_y_0_yyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 319);

            auto g_y_0_yyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 320);

            auto g_y_0_yyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 321);

            auto g_y_0_yyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 322);

            auto g_y_0_yyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 323);

            auto g_y_0_yyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 324);

            auto g_y_0_yyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 325);

            auto g_y_0_yyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 326);

            auto g_y_0_yyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 327);

            auto g_y_0_yyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 328);

            auto g_y_0_yyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 329);

            auto g_y_0_yyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 330);

            auto g_y_0_yyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 331);

            auto g_y_0_yyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 332);

            auto g_y_0_yyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 333);

            auto g_y_0_yyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 334);

            auto g_y_0_yyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 335);

            auto g_y_0_yyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 336);

            auto g_y_0_yyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 337);

            auto g_y_0_yyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 338);

            auto g_y_0_yyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 339);

            auto g_y_0_yyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 340);

            auto g_y_0_yyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 341);

            auto g_y_0_yyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 342);

            auto g_y_0_yyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 343);

            auto g_y_0_yyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 344);

            auto g_y_0_yyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 345);

            auto g_y_0_yyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 346);

            auto g_y_0_yyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 347);

            auto g_y_0_yyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 348);

            auto g_y_0_yyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 349);

            auto g_y_0_yyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 350);

            auto g_y_0_yyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 351);

            auto g_y_0_yyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 352);

            auto g_y_0_yyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 353);

            auto g_y_0_yyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 354);

            auto g_y_0_yyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 355);

            auto g_y_0_yyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 356);

            auto g_y_0_yyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 357);

            auto g_y_0_yyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 358);

            auto g_y_0_yyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 359);

            auto g_y_0_yyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 360);

            auto g_y_0_yyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 361);

            auto g_y_0_yyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 362);

            auto g_y_0_yyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 363);

            auto g_y_0_yzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 364);

            auto g_y_0_yzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 365);

            auto g_y_0_yzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 366);

            auto g_y_0_yzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 367);

            auto g_y_0_yzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 368);

            auto g_y_0_yzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 369);

            auto g_y_0_yzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 370);

            auto g_y_0_yzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 371);

            auto g_y_0_yzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 372);

            auto g_y_0_yzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 373);

            auto g_y_0_yzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 374);

            auto g_y_0_yzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 375);

            auto g_y_0_yzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 376);

            auto g_y_0_yzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 377);

            auto g_y_0_yzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 378);

            auto g_y_0_yzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 379);

            auto g_y_0_yzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 380);

            auto g_y_0_yzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 381);

            auto g_y_0_yzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 382);

            auto g_y_0_yzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 383);

            auto g_y_0_yzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 384);

            auto g_y_0_yzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 385);

            auto g_y_0_yzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 386);

            auto g_y_0_yzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 387);

            auto g_y_0_yzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 388);

            auto g_y_0_yzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 389);

            auto g_y_0_yzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 390);

            auto g_y_0_yzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 391);

            auto g_y_0_zzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 392);

            auto g_y_0_zzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 393);

            auto g_y_0_zzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 394);

            auto g_y_0_zzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 395);

            auto g_y_0_zzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 396);

            auto g_y_0_zzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 397);

            auto g_y_0_zzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 398);

            auto g_y_0_zzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 399);

            auto g_y_0_zzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 400);

            auto g_y_0_zzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 401);

            auto g_y_0_zzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 402);

            auto g_y_0_zzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 403);

            auto g_y_0_zzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 404);

            auto g_y_0_zzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 405);

            auto g_y_0_zzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 406);

            auto g_y_0_zzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 407);

            auto g_y_0_zzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 408);

            auto g_y_0_zzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 409);

            auto g_y_0_zzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 410);

            auto g_y_0_zzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 411);

            auto g_y_0_zzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 412);

            auto g_y_0_zzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 413);

            auto g_y_0_zzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 414);

            auto g_y_0_zzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 415);

            auto g_y_0_zzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 416);

            auto g_y_0_zzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 417);

            auto g_y_0_zzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 418);

            auto g_y_0_zzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 419);

            auto g_z_0_xxxx_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 0);

            auto g_z_0_xxxx_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 1);

            auto g_z_0_xxxx_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 2);

            auto g_z_0_xxxx_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 3);

            auto g_z_0_xxxx_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 4);

            auto g_z_0_xxxx_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 5);

            auto g_z_0_xxxx_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 6);

            auto g_z_0_xxxx_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 7);

            auto g_z_0_xxxx_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 8);

            auto g_z_0_xxxx_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 9);

            auto g_z_0_xxxx_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 10);

            auto g_z_0_xxxx_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 11);

            auto g_z_0_xxxx_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 12);

            auto g_z_0_xxxx_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 13);

            auto g_z_0_xxxx_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 14);

            auto g_z_0_xxxx_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 15);

            auto g_z_0_xxxx_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 16);

            auto g_z_0_xxxx_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 17);

            auto g_z_0_xxxx_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 18);

            auto g_z_0_xxxx_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 19);

            auto g_z_0_xxxx_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 20);

            auto g_z_0_xxxx_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 21);

            auto g_z_0_xxxx_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 22);

            auto g_z_0_xxxx_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 23);

            auto g_z_0_xxxx_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 24);

            auto g_z_0_xxxx_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 25);

            auto g_z_0_xxxx_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 26);

            auto g_z_0_xxxx_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 27);

            auto g_z_0_xxxy_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 28);

            auto g_z_0_xxxy_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 29);

            auto g_z_0_xxxy_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 30);

            auto g_z_0_xxxy_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 31);

            auto g_z_0_xxxy_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 32);

            auto g_z_0_xxxy_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 33);

            auto g_z_0_xxxy_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 34);

            auto g_z_0_xxxy_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 35);

            auto g_z_0_xxxy_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 36);

            auto g_z_0_xxxy_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 37);

            auto g_z_0_xxxy_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 38);

            auto g_z_0_xxxy_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 39);

            auto g_z_0_xxxy_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 40);

            auto g_z_0_xxxy_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 41);

            auto g_z_0_xxxy_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 42);

            auto g_z_0_xxxy_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 43);

            auto g_z_0_xxxy_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 44);

            auto g_z_0_xxxy_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 45);

            auto g_z_0_xxxy_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 46);

            auto g_z_0_xxxy_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 47);

            auto g_z_0_xxxy_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 48);

            auto g_z_0_xxxy_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 49);

            auto g_z_0_xxxy_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 50);

            auto g_z_0_xxxy_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 51);

            auto g_z_0_xxxy_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 52);

            auto g_z_0_xxxy_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 53);

            auto g_z_0_xxxy_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 54);

            auto g_z_0_xxxy_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 55);

            auto g_z_0_xxxz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 56);

            auto g_z_0_xxxz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 57);

            auto g_z_0_xxxz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 58);

            auto g_z_0_xxxz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 59);

            auto g_z_0_xxxz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 60);

            auto g_z_0_xxxz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 61);

            auto g_z_0_xxxz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 62);

            auto g_z_0_xxxz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 63);

            auto g_z_0_xxxz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 64);

            auto g_z_0_xxxz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 65);

            auto g_z_0_xxxz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 66);

            auto g_z_0_xxxz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 67);

            auto g_z_0_xxxz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 68);

            auto g_z_0_xxxz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 69);

            auto g_z_0_xxxz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 70);

            auto g_z_0_xxxz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 71);

            auto g_z_0_xxxz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 72);

            auto g_z_0_xxxz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 73);

            auto g_z_0_xxxz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 74);

            auto g_z_0_xxxz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 75);

            auto g_z_0_xxxz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 76);

            auto g_z_0_xxxz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 77);

            auto g_z_0_xxxz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 78);

            auto g_z_0_xxxz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 79);

            auto g_z_0_xxxz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 80);

            auto g_z_0_xxxz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 81);

            auto g_z_0_xxxz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 82);

            auto g_z_0_xxxz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 83);

            auto g_z_0_xxyy_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 84);

            auto g_z_0_xxyy_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 85);

            auto g_z_0_xxyy_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 86);

            auto g_z_0_xxyy_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 87);

            auto g_z_0_xxyy_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 88);

            auto g_z_0_xxyy_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 89);

            auto g_z_0_xxyy_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 90);

            auto g_z_0_xxyy_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 91);

            auto g_z_0_xxyy_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 92);

            auto g_z_0_xxyy_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 93);

            auto g_z_0_xxyy_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 94);

            auto g_z_0_xxyy_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 95);

            auto g_z_0_xxyy_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 96);

            auto g_z_0_xxyy_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 97);

            auto g_z_0_xxyy_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 98);

            auto g_z_0_xxyy_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 99);

            auto g_z_0_xxyy_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 100);

            auto g_z_0_xxyy_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 101);

            auto g_z_0_xxyy_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 102);

            auto g_z_0_xxyy_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 103);

            auto g_z_0_xxyy_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 104);

            auto g_z_0_xxyy_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 105);

            auto g_z_0_xxyy_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 106);

            auto g_z_0_xxyy_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 107);

            auto g_z_0_xxyy_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 108);

            auto g_z_0_xxyy_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 109);

            auto g_z_0_xxyy_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 110);

            auto g_z_0_xxyy_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 111);

            auto g_z_0_xxyz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 112);

            auto g_z_0_xxyz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 113);

            auto g_z_0_xxyz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 114);

            auto g_z_0_xxyz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 115);

            auto g_z_0_xxyz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 116);

            auto g_z_0_xxyz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 117);

            auto g_z_0_xxyz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 118);

            auto g_z_0_xxyz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 119);

            auto g_z_0_xxyz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 120);

            auto g_z_0_xxyz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 121);

            auto g_z_0_xxyz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 122);

            auto g_z_0_xxyz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 123);

            auto g_z_0_xxyz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 124);

            auto g_z_0_xxyz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 125);

            auto g_z_0_xxyz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 126);

            auto g_z_0_xxyz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 127);

            auto g_z_0_xxyz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 128);

            auto g_z_0_xxyz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 129);

            auto g_z_0_xxyz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 130);

            auto g_z_0_xxyz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 131);

            auto g_z_0_xxyz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 132);

            auto g_z_0_xxyz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 133);

            auto g_z_0_xxyz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 134);

            auto g_z_0_xxyz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 135);

            auto g_z_0_xxyz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 136);

            auto g_z_0_xxyz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 137);

            auto g_z_0_xxyz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 138);

            auto g_z_0_xxyz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 139);

            auto g_z_0_xxzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 140);

            auto g_z_0_xxzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 141);

            auto g_z_0_xxzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 142);

            auto g_z_0_xxzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 143);

            auto g_z_0_xxzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 144);

            auto g_z_0_xxzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 145);

            auto g_z_0_xxzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 146);

            auto g_z_0_xxzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 147);

            auto g_z_0_xxzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 148);

            auto g_z_0_xxzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 149);

            auto g_z_0_xxzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 150);

            auto g_z_0_xxzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 151);

            auto g_z_0_xxzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 152);

            auto g_z_0_xxzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 153);

            auto g_z_0_xxzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 154);

            auto g_z_0_xxzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 155);

            auto g_z_0_xxzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 156);

            auto g_z_0_xxzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 157);

            auto g_z_0_xxzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 158);

            auto g_z_0_xxzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 159);

            auto g_z_0_xxzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 160);

            auto g_z_0_xxzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 161);

            auto g_z_0_xxzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 162);

            auto g_z_0_xxzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 163);

            auto g_z_0_xxzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 164);

            auto g_z_0_xxzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 165);

            auto g_z_0_xxzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 166);

            auto g_z_0_xxzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 167);

            auto g_z_0_xyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 168);

            auto g_z_0_xyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 169);

            auto g_z_0_xyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 170);

            auto g_z_0_xyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 171);

            auto g_z_0_xyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 172);

            auto g_z_0_xyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 173);

            auto g_z_0_xyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 174);

            auto g_z_0_xyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 175);

            auto g_z_0_xyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 176);

            auto g_z_0_xyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 177);

            auto g_z_0_xyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 178);

            auto g_z_0_xyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 179);

            auto g_z_0_xyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 180);

            auto g_z_0_xyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 181);

            auto g_z_0_xyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 182);

            auto g_z_0_xyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 183);

            auto g_z_0_xyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 184);

            auto g_z_0_xyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 185);

            auto g_z_0_xyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 186);

            auto g_z_0_xyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 187);

            auto g_z_0_xyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 188);

            auto g_z_0_xyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 189);

            auto g_z_0_xyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 190);

            auto g_z_0_xyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 191);

            auto g_z_0_xyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 192);

            auto g_z_0_xyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 193);

            auto g_z_0_xyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 194);

            auto g_z_0_xyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 195);

            auto g_z_0_xyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 196);

            auto g_z_0_xyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 197);

            auto g_z_0_xyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 198);

            auto g_z_0_xyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 199);

            auto g_z_0_xyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 200);

            auto g_z_0_xyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 201);

            auto g_z_0_xyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 202);

            auto g_z_0_xyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 203);

            auto g_z_0_xyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 204);

            auto g_z_0_xyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 205);

            auto g_z_0_xyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 206);

            auto g_z_0_xyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 207);

            auto g_z_0_xyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 208);

            auto g_z_0_xyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 209);

            auto g_z_0_xyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 210);

            auto g_z_0_xyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 211);

            auto g_z_0_xyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 212);

            auto g_z_0_xyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 213);

            auto g_z_0_xyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 214);

            auto g_z_0_xyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 215);

            auto g_z_0_xyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 216);

            auto g_z_0_xyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 217);

            auto g_z_0_xyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 218);

            auto g_z_0_xyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 219);

            auto g_z_0_xyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 220);

            auto g_z_0_xyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 221);

            auto g_z_0_xyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 222);

            auto g_z_0_xyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 223);

            auto g_z_0_xyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 224);

            auto g_z_0_xyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 225);

            auto g_z_0_xyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 226);

            auto g_z_0_xyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 227);

            auto g_z_0_xyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 228);

            auto g_z_0_xyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 229);

            auto g_z_0_xyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 230);

            auto g_z_0_xyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 231);

            auto g_z_0_xyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 232);

            auto g_z_0_xyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 233);

            auto g_z_0_xyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 234);

            auto g_z_0_xyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 235);

            auto g_z_0_xyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 236);

            auto g_z_0_xyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 237);

            auto g_z_0_xyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 238);

            auto g_z_0_xyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 239);

            auto g_z_0_xyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 240);

            auto g_z_0_xyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 241);

            auto g_z_0_xyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 242);

            auto g_z_0_xyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 243);

            auto g_z_0_xyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 244);

            auto g_z_0_xyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 245);

            auto g_z_0_xyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 246);

            auto g_z_0_xyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 247);

            auto g_z_0_xyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 248);

            auto g_z_0_xyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 249);

            auto g_z_0_xyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 250);

            auto g_z_0_xyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 251);

            auto g_z_0_xzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 252);

            auto g_z_0_xzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 253);

            auto g_z_0_xzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 254);

            auto g_z_0_xzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 255);

            auto g_z_0_xzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 256);

            auto g_z_0_xzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 257);

            auto g_z_0_xzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 258);

            auto g_z_0_xzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 259);

            auto g_z_0_xzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 260);

            auto g_z_0_xzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 261);

            auto g_z_0_xzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 262);

            auto g_z_0_xzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 263);

            auto g_z_0_xzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 264);

            auto g_z_0_xzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 265);

            auto g_z_0_xzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 266);

            auto g_z_0_xzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 267);

            auto g_z_0_xzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 268);

            auto g_z_0_xzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 269);

            auto g_z_0_xzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 270);

            auto g_z_0_xzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 271);

            auto g_z_0_xzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 272);

            auto g_z_0_xzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 273);

            auto g_z_0_xzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 274);

            auto g_z_0_xzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 275);

            auto g_z_0_xzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 276);

            auto g_z_0_xzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 277);

            auto g_z_0_xzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 278);

            auto g_z_0_xzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 279);

            auto g_z_0_yyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 280);

            auto g_z_0_yyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 281);

            auto g_z_0_yyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 282);

            auto g_z_0_yyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 283);

            auto g_z_0_yyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 284);

            auto g_z_0_yyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 285);

            auto g_z_0_yyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 286);

            auto g_z_0_yyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 287);

            auto g_z_0_yyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 288);

            auto g_z_0_yyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 289);

            auto g_z_0_yyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 290);

            auto g_z_0_yyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 291);

            auto g_z_0_yyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 292);

            auto g_z_0_yyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 293);

            auto g_z_0_yyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 294);

            auto g_z_0_yyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 295);

            auto g_z_0_yyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 296);

            auto g_z_0_yyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 297);

            auto g_z_0_yyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 298);

            auto g_z_0_yyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 299);

            auto g_z_0_yyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 300);

            auto g_z_0_yyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 301);

            auto g_z_0_yyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 302);

            auto g_z_0_yyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 303);

            auto g_z_0_yyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 304);

            auto g_z_0_yyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 305);

            auto g_z_0_yyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 306);

            auto g_z_0_yyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 307);

            auto g_z_0_yyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 308);

            auto g_z_0_yyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 309);

            auto g_z_0_yyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 310);

            auto g_z_0_yyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 311);

            auto g_z_0_yyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 312);

            auto g_z_0_yyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 313);

            auto g_z_0_yyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 314);

            auto g_z_0_yyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 315);

            auto g_z_0_yyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 316);

            auto g_z_0_yyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 317);

            auto g_z_0_yyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 318);

            auto g_z_0_yyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 319);

            auto g_z_0_yyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 320);

            auto g_z_0_yyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 321);

            auto g_z_0_yyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 322);

            auto g_z_0_yyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 323);

            auto g_z_0_yyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 324);

            auto g_z_0_yyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 325);

            auto g_z_0_yyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 326);

            auto g_z_0_yyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 327);

            auto g_z_0_yyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 328);

            auto g_z_0_yyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 329);

            auto g_z_0_yyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 330);

            auto g_z_0_yyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 331);

            auto g_z_0_yyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 332);

            auto g_z_0_yyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 333);

            auto g_z_0_yyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 334);

            auto g_z_0_yyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 335);

            auto g_z_0_yyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 336);

            auto g_z_0_yyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 337);

            auto g_z_0_yyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 338);

            auto g_z_0_yyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 339);

            auto g_z_0_yyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 340);

            auto g_z_0_yyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 341);

            auto g_z_0_yyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 342);

            auto g_z_0_yyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 343);

            auto g_z_0_yyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 344);

            auto g_z_0_yyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 345);

            auto g_z_0_yyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 346);

            auto g_z_0_yyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 347);

            auto g_z_0_yyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 348);

            auto g_z_0_yyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 349);

            auto g_z_0_yyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 350);

            auto g_z_0_yyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 351);

            auto g_z_0_yyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 352);

            auto g_z_0_yyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 353);

            auto g_z_0_yyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 354);

            auto g_z_0_yyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 355);

            auto g_z_0_yyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 356);

            auto g_z_0_yyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 357);

            auto g_z_0_yyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 358);

            auto g_z_0_yyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 359);

            auto g_z_0_yyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 360);

            auto g_z_0_yyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 361);

            auto g_z_0_yyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 362);

            auto g_z_0_yyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 363);

            auto g_z_0_yzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 364);

            auto g_z_0_yzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 365);

            auto g_z_0_yzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 366);

            auto g_z_0_yzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 367);

            auto g_z_0_yzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 368);

            auto g_z_0_yzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 369);

            auto g_z_0_yzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 370);

            auto g_z_0_yzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 371);

            auto g_z_0_yzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 372);

            auto g_z_0_yzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 373);

            auto g_z_0_yzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 374);

            auto g_z_0_yzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 375);

            auto g_z_0_yzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 376);

            auto g_z_0_yzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 377);

            auto g_z_0_yzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 378);

            auto g_z_0_yzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 379);

            auto g_z_0_yzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 380);

            auto g_z_0_yzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 381);

            auto g_z_0_yzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 382);

            auto g_z_0_yzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 383);

            auto g_z_0_yzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 384);

            auto g_z_0_yzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 385);

            auto g_z_0_yzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 386);

            auto g_z_0_yzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 387);

            auto g_z_0_yzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 388);

            auto g_z_0_yzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 389);

            auto g_z_0_yzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 390);

            auto g_z_0_yzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 391);

            auto g_z_0_zzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 392);

            auto g_z_0_zzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 393);

            auto g_z_0_zzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 394);

            auto g_z_0_zzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 395);

            auto g_z_0_zzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 396);

            auto g_z_0_zzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 397);

            auto g_z_0_zzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 398);

            auto g_z_0_zzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 399);

            auto g_z_0_zzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 400);

            auto g_z_0_zzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 401);

            auto g_z_0_zzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 402);

            auto g_z_0_zzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 403);

            auto g_z_0_zzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 404);

            auto g_z_0_zzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 405);

            auto g_z_0_zzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 406);

            auto g_z_0_zzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 407);

            auto g_z_0_zzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 408);

            auto g_z_0_zzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 409);

            auto g_z_0_zzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 410);

            auto g_z_0_zzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 411);

            auto g_z_0_zzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 412);

            auto g_z_0_zzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 413);

            auto g_z_0_zzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 414);

            auto g_z_0_zzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 415);

            auto g_z_0_zzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 416);

            auto g_z_0_zzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 417);

            auto g_z_0_zzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 418);

            auto g_z_0_zzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 419);

            /// Set up components of auxilary buffer : SSGK

            const auto gk_geom_10_off = idx_geom_10_xxgk + (i * bcomps + j) * 540;

            auto g_x_0_xxxx_xxxxxxx = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxx_xxxxxxy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxx_xxxxxxz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxx_xxxxxyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxx_xxxxxyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxx_xxxxxzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxx_xxxxyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxx_xxxxyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxx_xxxxyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxx_xxxxzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxx_xxxyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxx_xxxyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxx_xxxyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxx_xxxyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxx_xxxzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxx_xxyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxx_xxyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxx_xxyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxx_xxyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxx_xxyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxx_xxzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxx_xyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxx_xyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxx_xyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxx_xyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxx_xyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxx_xyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxx_xzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxx_yyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxx_yyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxxx_yyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxx_yyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxx_yyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxx_yyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxx_yzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxx_zzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxy_xxxxxxx = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxy_xxxxxxy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxy_xxxxxxz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxy_xxxxxyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxy_xxxxxyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxy_xxxxxzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxxy_xxxxyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxy_xxxxyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxy_xxxxyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxxy_xxxxzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxy_xxxyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxy_xxxyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxxy_xxxyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxy_xxxyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxxy_xxxzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxy_xxyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxy_xxyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxy_xxyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxy_xxyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxy_xxyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxxy_xxzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxy_xyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxy_xyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxy_xyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xxxy_xyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxxy_xyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxxy_xyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxxy_xzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxxy_yyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxxy_yyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxxy_yyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxxy_yyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxxy_yyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxxy_yyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxxy_yzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxxy_zzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxxz_xxxxxxx = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxxz_xxxxxxy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxxz_xxxxxxz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxxz_xxxxxyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxxz_xxxxxyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxxz_xxxxxzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxxz_xxxxyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxxz_xxxxyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxxz_xxxxyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxxz_xxxxzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxxz_xxxyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxxz_xxxyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xxxz_xxxyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxxz_xxxyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxxz_xxxzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxxz_xxyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxxz_xxyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxxz_xxyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xxxz_xxyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxxz_xxyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxxz_xxzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxxz_xyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxxz_xyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxxz_xyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxxz_xyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxxz_xyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxxz_xyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxxz_xzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xxxz_yyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xxxz_yyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xxxz_yyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xxxz_yyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xxxz_yyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xxxz_yyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xxxz_yzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xxxz_zzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xxyy_xxxxxxx = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xxyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xxyy_xxxxxxz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xxyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xxyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xxyy_xxxxxzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xxyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xxyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xxyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xxyy_xxxxzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xxyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xxyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xxyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xxyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xxyy_xxxzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xxyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xxyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xxyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xxyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xxyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xxyy_xxzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xxyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xxyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xxyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xxyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xxyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xxyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xxyy_xzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xxyy_yyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xxyy_yyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xxyy_yyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xxyy_yyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_xxyy_yyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xxyy_yyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xxyy_yzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xxyy_zzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xxyz_xxxxxxx = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xxyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xxyz_xxxxxxz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_xxyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xxyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xxyz_xxxxxzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_xxyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xxyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xxyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_xxyz_xxxxzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xxyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xxyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xxyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_xxyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xxyz_xxxzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xxyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_xxyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xxyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_xxyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xxyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xxyz_xxzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_xxyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xxyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xxyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_x_0_xxyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_xxyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_xxyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_xxyz_xzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_xxyz_yyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_xxyz_yyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_xxyz_yyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_xxyz_yyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_xxyz_yyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_xxyz_yyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_xxyz_yzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_xxyz_zzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_xxzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_xxzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_xxzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_xxzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_xxzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_xxzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_xxzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_xxzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_xxzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_xxzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_xxzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_xxzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_xxzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_xxzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_xxzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_xxzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_xxzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_xxzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_xxzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_xxzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_xxzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_xxzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_xxzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_xxzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_xxzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_xxzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_xxzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_xxzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_xxzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_xxzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 209);

            auto g_x_0_xxzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_xxzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_xxzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_xxzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_xxzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_xxzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_xyyy_xxxxxxx = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_xyyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_xyyy_xxxxxxz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_xyyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_xyyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_xyyy_xxxxxzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_xyyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_xyyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 223);

            auto g_x_0_xyyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 224);

            auto g_x_0_xyyy_xxxxzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_xyyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_xyyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_xyyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_xyyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 229);

            auto g_x_0_xyyy_xxxzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 230);

            auto g_x_0_xyyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_xyyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_xyyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 233);

            auto g_x_0_xyyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_xyyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_xyyy_xxzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_xyyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_xyyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_xyyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 239);

            auto g_x_0_xyyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_xyyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_xyyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_xyyy_xzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_xyyy_yyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_xyyy_yyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_xyyy_yyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_xyyy_yyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_xyyy_yyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_xyyy_yyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_xyyy_yzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_xyyy_zzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 251);

            auto g_x_0_xyyz_xxxxxxx = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 252);

            auto g_x_0_xyyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_xyyz_xxxxxxz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 254);

            auto g_x_0_xyyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_xyyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_xyyz_xxxxxzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 257);

            auto g_x_0_xyyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_xyyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 259);

            auto g_x_0_xyyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_xyyz_xxxxzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 261);

            auto g_x_0_xyyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_xyyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_xyyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_xyyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_xyyz_xxxzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 266);

            auto g_x_0_xyyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_xyyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_xyyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 269);

            auto g_x_0_xyyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_xyyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_xyyz_xxzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 272);

            auto g_x_0_xyyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 273);

            auto g_x_0_xyyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_xyyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 275);

            auto g_x_0_xyyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_xyyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_xyyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 278);

            auto g_x_0_xyyz_xzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 279);

            auto g_x_0_xyyz_yyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 280);

            auto g_x_0_xyyz_yyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 281);

            auto g_x_0_xyyz_yyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 282);

            auto g_x_0_xyyz_yyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 283);

            auto g_x_0_xyyz_yyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 284);

            auto g_x_0_xyyz_yyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 285);

            auto g_x_0_xyyz_yzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 286);

            auto g_x_0_xyyz_zzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 287);

            auto g_x_0_xyzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 288);

            auto g_x_0_xyzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 289);

            auto g_x_0_xyzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 290);

            auto g_x_0_xyzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 291);

            auto g_x_0_xyzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 292);

            auto g_x_0_xyzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 293);

            auto g_x_0_xyzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 294);

            auto g_x_0_xyzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 295);

            auto g_x_0_xyzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 296);

            auto g_x_0_xyzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 297);

            auto g_x_0_xyzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 298);

            auto g_x_0_xyzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 299);

            auto g_x_0_xyzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 300);

            auto g_x_0_xyzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 301);

            auto g_x_0_xyzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 302);

            auto g_x_0_xyzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 303);

            auto g_x_0_xyzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 304);

            auto g_x_0_xyzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 305);

            auto g_x_0_xyzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 306);

            auto g_x_0_xyzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 307);

            auto g_x_0_xyzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 308);

            auto g_x_0_xyzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 309);

            auto g_x_0_xyzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 310);

            auto g_x_0_xyzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 311);

            auto g_x_0_xyzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 312);

            auto g_x_0_xyzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 313);

            auto g_x_0_xyzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 314);

            auto g_x_0_xyzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 315);

            auto g_x_0_xyzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 316);

            auto g_x_0_xyzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 317);

            auto g_x_0_xyzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 318);

            auto g_x_0_xyzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 319);

            auto g_x_0_xyzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 320);

            auto g_x_0_xyzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 321);

            auto g_x_0_xyzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 322);

            auto g_x_0_xyzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 323);

            auto g_x_0_xzzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 324);

            auto g_x_0_xzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 325);

            auto g_x_0_xzzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 326);

            auto g_x_0_xzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 327);

            auto g_x_0_xzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 328);

            auto g_x_0_xzzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 329);

            auto g_x_0_xzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 330);

            auto g_x_0_xzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 331);

            auto g_x_0_xzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 332);

            auto g_x_0_xzzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 333);

            auto g_x_0_xzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 334);

            auto g_x_0_xzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 335);

            auto g_x_0_xzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 336);

            auto g_x_0_xzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 337);

            auto g_x_0_xzzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 338);

            auto g_x_0_xzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 339);

            auto g_x_0_xzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 340);

            auto g_x_0_xzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 341);

            auto g_x_0_xzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 342);

            auto g_x_0_xzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 343);

            auto g_x_0_xzzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 344);

            auto g_x_0_xzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 345);

            auto g_x_0_xzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 346);

            auto g_x_0_xzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 347);

            auto g_x_0_xzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 348);

            auto g_x_0_xzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 349);

            auto g_x_0_xzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 350);

            auto g_x_0_xzzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 351);

            auto g_x_0_xzzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 352);

            auto g_x_0_xzzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 353);

            auto g_x_0_xzzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 354);

            auto g_x_0_xzzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 355);

            auto g_x_0_xzzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 356);

            auto g_x_0_xzzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 357);

            auto g_x_0_xzzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 358);

            auto g_x_0_xzzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 359);

            auto g_x_0_yyyy_xxxxxxx = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 360);

            auto g_x_0_yyyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 361);

            auto g_x_0_yyyy_xxxxxxz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 362);

            auto g_x_0_yyyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 363);

            auto g_x_0_yyyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 364);

            auto g_x_0_yyyy_xxxxxzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 365);

            auto g_x_0_yyyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 366);

            auto g_x_0_yyyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 367);

            auto g_x_0_yyyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 368);

            auto g_x_0_yyyy_xxxxzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 369);

            auto g_x_0_yyyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 370);

            auto g_x_0_yyyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 371);

            auto g_x_0_yyyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 372);

            auto g_x_0_yyyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 373);

            auto g_x_0_yyyy_xxxzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 374);

            auto g_x_0_yyyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 375);

            auto g_x_0_yyyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 376);

            auto g_x_0_yyyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 377);

            auto g_x_0_yyyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 378);

            auto g_x_0_yyyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 379);

            auto g_x_0_yyyy_xxzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 380);

            auto g_x_0_yyyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 381);

            auto g_x_0_yyyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 382);

            auto g_x_0_yyyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 383);

            auto g_x_0_yyyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 384);

            auto g_x_0_yyyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 385);

            auto g_x_0_yyyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 386);

            auto g_x_0_yyyy_xzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 387);

            auto g_x_0_yyyy_yyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 388);

            auto g_x_0_yyyy_yyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 389);

            auto g_x_0_yyyy_yyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 390);

            auto g_x_0_yyyy_yyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 391);

            auto g_x_0_yyyy_yyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 392);

            auto g_x_0_yyyy_yyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 393);

            auto g_x_0_yyyy_yzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 394);

            auto g_x_0_yyyy_zzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 395);

            auto g_x_0_yyyz_xxxxxxx = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 396);

            auto g_x_0_yyyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 397);

            auto g_x_0_yyyz_xxxxxxz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 398);

            auto g_x_0_yyyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 399);

            auto g_x_0_yyyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 400);

            auto g_x_0_yyyz_xxxxxzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 401);

            auto g_x_0_yyyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 402);

            auto g_x_0_yyyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 403);

            auto g_x_0_yyyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 404);

            auto g_x_0_yyyz_xxxxzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 405);

            auto g_x_0_yyyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 406);

            auto g_x_0_yyyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 407);

            auto g_x_0_yyyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 408);

            auto g_x_0_yyyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 409);

            auto g_x_0_yyyz_xxxzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 410);

            auto g_x_0_yyyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 411);

            auto g_x_0_yyyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 412);

            auto g_x_0_yyyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 413);

            auto g_x_0_yyyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 414);

            auto g_x_0_yyyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 415);

            auto g_x_0_yyyz_xxzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 416);

            auto g_x_0_yyyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 417);

            auto g_x_0_yyyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 418);

            auto g_x_0_yyyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 419);

            auto g_x_0_yyyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 420);

            auto g_x_0_yyyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 421);

            auto g_x_0_yyyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 422);

            auto g_x_0_yyyz_xzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 423);

            auto g_x_0_yyyz_yyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 424);

            auto g_x_0_yyyz_yyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 425);

            auto g_x_0_yyyz_yyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 426);

            auto g_x_0_yyyz_yyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 427);

            auto g_x_0_yyyz_yyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 428);

            auto g_x_0_yyyz_yyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 429);

            auto g_x_0_yyyz_yzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 430);

            auto g_x_0_yyyz_zzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 431);

            auto g_x_0_yyzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 432);

            auto g_x_0_yyzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 433);

            auto g_x_0_yyzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 434);

            auto g_x_0_yyzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 435);

            auto g_x_0_yyzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 436);

            auto g_x_0_yyzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 437);

            auto g_x_0_yyzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 438);

            auto g_x_0_yyzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 439);

            auto g_x_0_yyzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 440);

            auto g_x_0_yyzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 441);

            auto g_x_0_yyzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 442);

            auto g_x_0_yyzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 443);

            auto g_x_0_yyzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 444);

            auto g_x_0_yyzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 445);

            auto g_x_0_yyzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 446);

            auto g_x_0_yyzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 447);

            auto g_x_0_yyzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 448);

            auto g_x_0_yyzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 449);

            auto g_x_0_yyzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 450);

            auto g_x_0_yyzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 451);

            auto g_x_0_yyzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 452);

            auto g_x_0_yyzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 453);

            auto g_x_0_yyzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 454);

            auto g_x_0_yyzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 455);

            auto g_x_0_yyzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 456);

            auto g_x_0_yyzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 457);

            auto g_x_0_yyzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 458);

            auto g_x_0_yyzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 459);

            auto g_x_0_yyzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 460);

            auto g_x_0_yyzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 461);

            auto g_x_0_yyzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 462);

            auto g_x_0_yyzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 463);

            auto g_x_0_yyzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 464);

            auto g_x_0_yyzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 465);

            auto g_x_0_yyzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 466);

            auto g_x_0_yyzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 467);

            auto g_x_0_yzzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 468);

            auto g_x_0_yzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 469);

            auto g_x_0_yzzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 470);

            auto g_x_0_yzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 471);

            auto g_x_0_yzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 472);

            auto g_x_0_yzzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 473);

            auto g_x_0_yzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 474);

            auto g_x_0_yzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 475);

            auto g_x_0_yzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 476);

            auto g_x_0_yzzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 477);

            auto g_x_0_yzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 478);

            auto g_x_0_yzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 479);

            auto g_x_0_yzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 480);

            auto g_x_0_yzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 481);

            auto g_x_0_yzzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 482);

            auto g_x_0_yzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 483);

            auto g_x_0_yzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 484);

            auto g_x_0_yzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 485);

            auto g_x_0_yzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 486);

            auto g_x_0_yzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 487);

            auto g_x_0_yzzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 488);

            auto g_x_0_yzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 489);

            auto g_x_0_yzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 490);

            auto g_x_0_yzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 491);

            auto g_x_0_yzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 492);

            auto g_x_0_yzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 493);

            auto g_x_0_yzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 494);

            auto g_x_0_yzzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 495);

            auto g_x_0_yzzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 496);

            auto g_x_0_yzzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 497);

            auto g_x_0_yzzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 498);

            auto g_x_0_yzzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 499);

            auto g_x_0_yzzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 500);

            auto g_x_0_yzzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 501);

            auto g_x_0_yzzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 502);

            auto g_x_0_yzzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 503);

            auto g_x_0_zzzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 504);

            auto g_x_0_zzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 505);

            auto g_x_0_zzzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 506);

            auto g_x_0_zzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 507);

            auto g_x_0_zzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 508);

            auto g_x_0_zzzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 509);

            auto g_x_0_zzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 510);

            auto g_x_0_zzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 511);

            auto g_x_0_zzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 512);

            auto g_x_0_zzzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 513);

            auto g_x_0_zzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 514);

            auto g_x_0_zzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 515);

            auto g_x_0_zzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 516);

            auto g_x_0_zzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 517);

            auto g_x_0_zzzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 518);

            auto g_x_0_zzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 519);

            auto g_x_0_zzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 520);

            auto g_x_0_zzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 521);

            auto g_x_0_zzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 522);

            auto g_x_0_zzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 523);

            auto g_x_0_zzzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 524);

            auto g_x_0_zzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 525);

            auto g_x_0_zzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 526);

            auto g_x_0_zzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 527);

            auto g_x_0_zzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 528);

            auto g_x_0_zzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 529);

            auto g_x_0_zzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 530);

            auto g_x_0_zzzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 531);

            auto g_x_0_zzzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 532);

            auto g_x_0_zzzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 533);

            auto g_x_0_zzzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 534);

            auto g_x_0_zzzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 535);

            auto g_x_0_zzzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 536);

            auto g_x_0_zzzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 537);

            auto g_x_0_zzzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 538);

            auto g_x_0_zzzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 0 * acomps * bcomps + 539);

            auto g_y_0_xxxx_xxxxxxx = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 0);

            auto g_y_0_xxxx_xxxxxxy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 1);

            auto g_y_0_xxxx_xxxxxxz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 2);

            auto g_y_0_xxxx_xxxxxyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 3);

            auto g_y_0_xxxx_xxxxxyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 4);

            auto g_y_0_xxxx_xxxxxzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 5);

            auto g_y_0_xxxx_xxxxyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 6);

            auto g_y_0_xxxx_xxxxyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 7);

            auto g_y_0_xxxx_xxxxyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 8);

            auto g_y_0_xxxx_xxxxzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 9);

            auto g_y_0_xxxx_xxxyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 10);

            auto g_y_0_xxxx_xxxyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 11);

            auto g_y_0_xxxx_xxxyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 12);

            auto g_y_0_xxxx_xxxyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 13);

            auto g_y_0_xxxx_xxxzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 14);

            auto g_y_0_xxxx_xxyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 15);

            auto g_y_0_xxxx_xxyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 16);

            auto g_y_0_xxxx_xxyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 17);

            auto g_y_0_xxxx_xxyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 18);

            auto g_y_0_xxxx_xxyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 19);

            auto g_y_0_xxxx_xxzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 20);

            auto g_y_0_xxxx_xyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 21);

            auto g_y_0_xxxx_xyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 22);

            auto g_y_0_xxxx_xyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 23);

            auto g_y_0_xxxx_xyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 24);

            auto g_y_0_xxxx_xyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 25);

            auto g_y_0_xxxx_xyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 26);

            auto g_y_0_xxxx_xzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 27);

            auto g_y_0_xxxx_yyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 28);

            auto g_y_0_xxxx_yyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 29);

            auto g_y_0_xxxx_yyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 30);

            auto g_y_0_xxxx_yyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 31);

            auto g_y_0_xxxx_yyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 32);

            auto g_y_0_xxxx_yyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 33);

            auto g_y_0_xxxx_yzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 34);

            auto g_y_0_xxxx_zzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 35);

            auto g_y_0_xxxy_xxxxxxx = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 36);

            auto g_y_0_xxxy_xxxxxxy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 37);

            auto g_y_0_xxxy_xxxxxxz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 38);

            auto g_y_0_xxxy_xxxxxyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 39);

            auto g_y_0_xxxy_xxxxxyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 40);

            auto g_y_0_xxxy_xxxxxzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 41);

            auto g_y_0_xxxy_xxxxyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 42);

            auto g_y_0_xxxy_xxxxyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 43);

            auto g_y_0_xxxy_xxxxyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 44);

            auto g_y_0_xxxy_xxxxzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 45);

            auto g_y_0_xxxy_xxxyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 46);

            auto g_y_0_xxxy_xxxyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 47);

            auto g_y_0_xxxy_xxxyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 48);

            auto g_y_0_xxxy_xxxyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 49);

            auto g_y_0_xxxy_xxxzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 50);

            auto g_y_0_xxxy_xxyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 51);

            auto g_y_0_xxxy_xxyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 52);

            auto g_y_0_xxxy_xxyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 53);

            auto g_y_0_xxxy_xxyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 54);

            auto g_y_0_xxxy_xxyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 55);

            auto g_y_0_xxxy_xxzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 56);

            auto g_y_0_xxxy_xyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 57);

            auto g_y_0_xxxy_xyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 58);

            auto g_y_0_xxxy_xyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 59);

            auto g_y_0_xxxy_xyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 60);

            auto g_y_0_xxxy_xyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 61);

            auto g_y_0_xxxy_xyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 62);

            auto g_y_0_xxxy_xzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 63);

            auto g_y_0_xxxy_yyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 64);

            auto g_y_0_xxxy_yyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 65);

            auto g_y_0_xxxy_yyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 66);

            auto g_y_0_xxxy_yyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 67);

            auto g_y_0_xxxy_yyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 68);

            auto g_y_0_xxxy_yyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 69);

            auto g_y_0_xxxy_yzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 70);

            auto g_y_0_xxxy_zzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 71);

            auto g_y_0_xxxz_xxxxxxx = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 72);

            auto g_y_0_xxxz_xxxxxxy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 73);

            auto g_y_0_xxxz_xxxxxxz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 74);

            auto g_y_0_xxxz_xxxxxyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 75);

            auto g_y_0_xxxz_xxxxxyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 76);

            auto g_y_0_xxxz_xxxxxzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 77);

            auto g_y_0_xxxz_xxxxyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 78);

            auto g_y_0_xxxz_xxxxyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 79);

            auto g_y_0_xxxz_xxxxyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 80);

            auto g_y_0_xxxz_xxxxzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 81);

            auto g_y_0_xxxz_xxxyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 82);

            auto g_y_0_xxxz_xxxyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 83);

            auto g_y_0_xxxz_xxxyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 84);

            auto g_y_0_xxxz_xxxyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 85);

            auto g_y_0_xxxz_xxxzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 86);

            auto g_y_0_xxxz_xxyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 87);

            auto g_y_0_xxxz_xxyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 88);

            auto g_y_0_xxxz_xxyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 89);

            auto g_y_0_xxxz_xxyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 90);

            auto g_y_0_xxxz_xxyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 91);

            auto g_y_0_xxxz_xxzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 92);

            auto g_y_0_xxxz_xyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 93);

            auto g_y_0_xxxz_xyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 94);

            auto g_y_0_xxxz_xyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 95);

            auto g_y_0_xxxz_xyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 96);

            auto g_y_0_xxxz_xyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 97);

            auto g_y_0_xxxz_xyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 98);

            auto g_y_0_xxxz_xzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 99);

            auto g_y_0_xxxz_yyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 100);

            auto g_y_0_xxxz_yyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 101);

            auto g_y_0_xxxz_yyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 102);

            auto g_y_0_xxxz_yyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 103);

            auto g_y_0_xxxz_yyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 104);

            auto g_y_0_xxxz_yyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 105);

            auto g_y_0_xxxz_yzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 106);

            auto g_y_0_xxxz_zzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 107);

            auto g_y_0_xxyy_xxxxxxx = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 108);

            auto g_y_0_xxyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 109);

            auto g_y_0_xxyy_xxxxxxz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 110);

            auto g_y_0_xxyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 111);

            auto g_y_0_xxyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 112);

            auto g_y_0_xxyy_xxxxxzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 113);

            auto g_y_0_xxyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 114);

            auto g_y_0_xxyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 115);

            auto g_y_0_xxyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 116);

            auto g_y_0_xxyy_xxxxzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 117);

            auto g_y_0_xxyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 118);

            auto g_y_0_xxyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 119);

            auto g_y_0_xxyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 120);

            auto g_y_0_xxyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 121);

            auto g_y_0_xxyy_xxxzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 122);

            auto g_y_0_xxyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 123);

            auto g_y_0_xxyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 124);

            auto g_y_0_xxyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 125);

            auto g_y_0_xxyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 126);

            auto g_y_0_xxyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 127);

            auto g_y_0_xxyy_xxzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 128);

            auto g_y_0_xxyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 129);

            auto g_y_0_xxyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 130);

            auto g_y_0_xxyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 131);

            auto g_y_0_xxyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 132);

            auto g_y_0_xxyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 133);

            auto g_y_0_xxyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 134);

            auto g_y_0_xxyy_xzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 135);

            auto g_y_0_xxyy_yyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 136);

            auto g_y_0_xxyy_yyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 137);

            auto g_y_0_xxyy_yyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 138);

            auto g_y_0_xxyy_yyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 139);

            auto g_y_0_xxyy_yyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 140);

            auto g_y_0_xxyy_yyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 141);

            auto g_y_0_xxyy_yzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 142);

            auto g_y_0_xxyy_zzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 143);

            auto g_y_0_xxyz_xxxxxxx = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 144);

            auto g_y_0_xxyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 145);

            auto g_y_0_xxyz_xxxxxxz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 146);

            auto g_y_0_xxyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 147);

            auto g_y_0_xxyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 148);

            auto g_y_0_xxyz_xxxxxzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 149);

            auto g_y_0_xxyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 150);

            auto g_y_0_xxyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 151);

            auto g_y_0_xxyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 152);

            auto g_y_0_xxyz_xxxxzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 153);

            auto g_y_0_xxyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 154);

            auto g_y_0_xxyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 155);

            auto g_y_0_xxyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 156);

            auto g_y_0_xxyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 157);

            auto g_y_0_xxyz_xxxzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 158);

            auto g_y_0_xxyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 159);

            auto g_y_0_xxyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 160);

            auto g_y_0_xxyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 161);

            auto g_y_0_xxyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 162);

            auto g_y_0_xxyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 163);

            auto g_y_0_xxyz_xxzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 164);

            auto g_y_0_xxyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 165);

            auto g_y_0_xxyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 166);

            auto g_y_0_xxyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 167);

            auto g_y_0_xxyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 168);

            auto g_y_0_xxyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 169);

            auto g_y_0_xxyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 170);

            auto g_y_0_xxyz_xzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 171);

            auto g_y_0_xxyz_yyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 172);

            auto g_y_0_xxyz_yyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 173);

            auto g_y_0_xxyz_yyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 174);

            auto g_y_0_xxyz_yyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 175);

            auto g_y_0_xxyz_yyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 176);

            auto g_y_0_xxyz_yyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 177);

            auto g_y_0_xxyz_yzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 178);

            auto g_y_0_xxyz_zzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 179);

            auto g_y_0_xxzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 180);

            auto g_y_0_xxzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 181);

            auto g_y_0_xxzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 182);

            auto g_y_0_xxzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 183);

            auto g_y_0_xxzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 184);

            auto g_y_0_xxzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 185);

            auto g_y_0_xxzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 186);

            auto g_y_0_xxzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 187);

            auto g_y_0_xxzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 188);

            auto g_y_0_xxzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 189);

            auto g_y_0_xxzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 190);

            auto g_y_0_xxzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 191);

            auto g_y_0_xxzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 192);

            auto g_y_0_xxzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 193);

            auto g_y_0_xxzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 194);

            auto g_y_0_xxzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 195);

            auto g_y_0_xxzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 196);

            auto g_y_0_xxzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 197);

            auto g_y_0_xxzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 198);

            auto g_y_0_xxzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 199);

            auto g_y_0_xxzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 200);

            auto g_y_0_xxzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 201);

            auto g_y_0_xxzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 202);

            auto g_y_0_xxzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 203);

            auto g_y_0_xxzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 204);

            auto g_y_0_xxzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 205);

            auto g_y_0_xxzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 206);

            auto g_y_0_xxzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 207);

            auto g_y_0_xxzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 208);

            auto g_y_0_xxzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 209);

            auto g_y_0_xxzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 210);

            auto g_y_0_xxzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 211);

            auto g_y_0_xxzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 212);

            auto g_y_0_xxzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 213);

            auto g_y_0_xxzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 214);

            auto g_y_0_xxzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 215);

            auto g_y_0_xyyy_xxxxxxx = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 216);

            auto g_y_0_xyyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 217);

            auto g_y_0_xyyy_xxxxxxz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 218);

            auto g_y_0_xyyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 219);

            auto g_y_0_xyyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 220);

            auto g_y_0_xyyy_xxxxxzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 221);

            auto g_y_0_xyyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 222);

            auto g_y_0_xyyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 223);

            auto g_y_0_xyyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 224);

            auto g_y_0_xyyy_xxxxzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 225);

            auto g_y_0_xyyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 226);

            auto g_y_0_xyyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 227);

            auto g_y_0_xyyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 228);

            auto g_y_0_xyyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 229);

            auto g_y_0_xyyy_xxxzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 230);

            auto g_y_0_xyyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 231);

            auto g_y_0_xyyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 232);

            auto g_y_0_xyyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 233);

            auto g_y_0_xyyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 234);

            auto g_y_0_xyyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 235);

            auto g_y_0_xyyy_xxzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 236);

            auto g_y_0_xyyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 237);

            auto g_y_0_xyyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 238);

            auto g_y_0_xyyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 239);

            auto g_y_0_xyyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 240);

            auto g_y_0_xyyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 241);

            auto g_y_0_xyyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 242);

            auto g_y_0_xyyy_xzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 243);

            auto g_y_0_xyyy_yyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 244);

            auto g_y_0_xyyy_yyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 245);

            auto g_y_0_xyyy_yyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 246);

            auto g_y_0_xyyy_yyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 247);

            auto g_y_0_xyyy_yyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 248);

            auto g_y_0_xyyy_yyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 249);

            auto g_y_0_xyyy_yzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 250);

            auto g_y_0_xyyy_zzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 251);

            auto g_y_0_xyyz_xxxxxxx = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 252);

            auto g_y_0_xyyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 253);

            auto g_y_0_xyyz_xxxxxxz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 254);

            auto g_y_0_xyyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 255);

            auto g_y_0_xyyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 256);

            auto g_y_0_xyyz_xxxxxzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 257);

            auto g_y_0_xyyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 258);

            auto g_y_0_xyyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 259);

            auto g_y_0_xyyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 260);

            auto g_y_0_xyyz_xxxxzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 261);

            auto g_y_0_xyyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 262);

            auto g_y_0_xyyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 263);

            auto g_y_0_xyyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 264);

            auto g_y_0_xyyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 265);

            auto g_y_0_xyyz_xxxzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 266);

            auto g_y_0_xyyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 267);

            auto g_y_0_xyyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 268);

            auto g_y_0_xyyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 269);

            auto g_y_0_xyyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 270);

            auto g_y_0_xyyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 271);

            auto g_y_0_xyyz_xxzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 272);

            auto g_y_0_xyyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 273);

            auto g_y_0_xyyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 274);

            auto g_y_0_xyyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 275);

            auto g_y_0_xyyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 276);

            auto g_y_0_xyyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 277);

            auto g_y_0_xyyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 278);

            auto g_y_0_xyyz_xzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 279);

            auto g_y_0_xyyz_yyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 280);

            auto g_y_0_xyyz_yyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 281);

            auto g_y_0_xyyz_yyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 282);

            auto g_y_0_xyyz_yyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 283);

            auto g_y_0_xyyz_yyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 284);

            auto g_y_0_xyyz_yyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 285);

            auto g_y_0_xyyz_yzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 286);

            auto g_y_0_xyyz_zzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 287);

            auto g_y_0_xyzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 288);

            auto g_y_0_xyzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 289);

            auto g_y_0_xyzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 290);

            auto g_y_0_xyzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 291);

            auto g_y_0_xyzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 292);

            auto g_y_0_xyzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 293);

            auto g_y_0_xyzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 294);

            auto g_y_0_xyzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 295);

            auto g_y_0_xyzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 296);

            auto g_y_0_xyzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 297);

            auto g_y_0_xyzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 298);

            auto g_y_0_xyzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 299);

            auto g_y_0_xyzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 300);

            auto g_y_0_xyzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 301);

            auto g_y_0_xyzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 302);

            auto g_y_0_xyzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 303);

            auto g_y_0_xyzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 304);

            auto g_y_0_xyzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 305);

            auto g_y_0_xyzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 306);

            auto g_y_0_xyzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 307);

            auto g_y_0_xyzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 308);

            auto g_y_0_xyzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 309);

            auto g_y_0_xyzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 310);

            auto g_y_0_xyzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 311);

            auto g_y_0_xyzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 312);

            auto g_y_0_xyzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 313);

            auto g_y_0_xyzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 314);

            auto g_y_0_xyzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 315);

            auto g_y_0_xyzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 316);

            auto g_y_0_xyzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 317);

            auto g_y_0_xyzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 318);

            auto g_y_0_xyzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 319);

            auto g_y_0_xyzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 320);

            auto g_y_0_xyzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 321);

            auto g_y_0_xyzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 322);

            auto g_y_0_xyzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 323);

            auto g_y_0_xzzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 324);

            auto g_y_0_xzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 325);

            auto g_y_0_xzzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 326);

            auto g_y_0_xzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 327);

            auto g_y_0_xzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 328);

            auto g_y_0_xzzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 329);

            auto g_y_0_xzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 330);

            auto g_y_0_xzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 331);

            auto g_y_0_xzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 332);

            auto g_y_0_xzzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 333);

            auto g_y_0_xzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 334);

            auto g_y_0_xzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 335);

            auto g_y_0_xzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 336);

            auto g_y_0_xzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 337);

            auto g_y_0_xzzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 338);

            auto g_y_0_xzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 339);

            auto g_y_0_xzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 340);

            auto g_y_0_xzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 341);

            auto g_y_0_xzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 342);

            auto g_y_0_xzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 343);

            auto g_y_0_xzzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 344);

            auto g_y_0_xzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 345);

            auto g_y_0_xzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 346);

            auto g_y_0_xzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 347);

            auto g_y_0_xzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 348);

            auto g_y_0_xzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 349);

            auto g_y_0_xzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 350);

            auto g_y_0_xzzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 351);

            auto g_y_0_xzzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 352);

            auto g_y_0_xzzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 353);

            auto g_y_0_xzzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 354);

            auto g_y_0_xzzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 355);

            auto g_y_0_xzzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 356);

            auto g_y_0_xzzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 357);

            auto g_y_0_xzzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 358);

            auto g_y_0_xzzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 359);

            auto g_y_0_yyyy_xxxxxxx = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 360);

            auto g_y_0_yyyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 361);

            auto g_y_0_yyyy_xxxxxxz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 362);

            auto g_y_0_yyyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 363);

            auto g_y_0_yyyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 364);

            auto g_y_0_yyyy_xxxxxzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 365);

            auto g_y_0_yyyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 366);

            auto g_y_0_yyyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 367);

            auto g_y_0_yyyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 368);

            auto g_y_0_yyyy_xxxxzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 369);

            auto g_y_0_yyyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 370);

            auto g_y_0_yyyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 371);

            auto g_y_0_yyyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 372);

            auto g_y_0_yyyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 373);

            auto g_y_0_yyyy_xxxzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 374);

            auto g_y_0_yyyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 375);

            auto g_y_0_yyyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 376);

            auto g_y_0_yyyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 377);

            auto g_y_0_yyyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 378);

            auto g_y_0_yyyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 379);

            auto g_y_0_yyyy_xxzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 380);

            auto g_y_0_yyyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 381);

            auto g_y_0_yyyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 382);

            auto g_y_0_yyyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 383);

            auto g_y_0_yyyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 384);

            auto g_y_0_yyyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 385);

            auto g_y_0_yyyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 386);

            auto g_y_0_yyyy_xzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 387);

            auto g_y_0_yyyy_yyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 388);

            auto g_y_0_yyyy_yyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 389);

            auto g_y_0_yyyy_yyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 390);

            auto g_y_0_yyyy_yyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 391);

            auto g_y_0_yyyy_yyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 392);

            auto g_y_0_yyyy_yyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 393);

            auto g_y_0_yyyy_yzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 394);

            auto g_y_0_yyyy_zzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 395);

            auto g_y_0_yyyz_xxxxxxx = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 396);

            auto g_y_0_yyyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 397);

            auto g_y_0_yyyz_xxxxxxz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 398);

            auto g_y_0_yyyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 399);

            auto g_y_0_yyyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 400);

            auto g_y_0_yyyz_xxxxxzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 401);

            auto g_y_0_yyyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 402);

            auto g_y_0_yyyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 403);

            auto g_y_0_yyyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 404);

            auto g_y_0_yyyz_xxxxzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 405);

            auto g_y_0_yyyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 406);

            auto g_y_0_yyyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 407);

            auto g_y_0_yyyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 408);

            auto g_y_0_yyyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 409);

            auto g_y_0_yyyz_xxxzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 410);

            auto g_y_0_yyyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 411);

            auto g_y_0_yyyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 412);

            auto g_y_0_yyyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 413);

            auto g_y_0_yyyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 414);

            auto g_y_0_yyyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 415);

            auto g_y_0_yyyz_xxzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 416);

            auto g_y_0_yyyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 417);

            auto g_y_0_yyyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 418);

            auto g_y_0_yyyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 419);

            auto g_y_0_yyyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 420);

            auto g_y_0_yyyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 421);

            auto g_y_0_yyyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 422);

            auto g_y_0_yyyz_xzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 423);

            auto g_y_0_yyyz_yyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 424);

            auto g_y_0_yyyz_yyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 425);

            auto g_y_0_yyyz_yyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 426);

            auto g_y_0_yyyz_yyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 427);

            auto g_y_0_yyyz_yyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 428);

            auto g_y_0_yyyz_yyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 429);

            auto g_y_0_yyyz_yzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 430);

            auto g_y_0_yyyz_zzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 431);

            auto g_y_0_yyzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 432);

            auto g_y_0_yyzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 433);

            auto g_y_0_yyzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 434);

            auto g_y_0_yyzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 435);

            auto g_y_0_yyzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 436);

            auto g_y_0_yyzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 437);

            auto g_y_0_yyzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 438);

            auto g_y_0_yyzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 439);

            auto g_y_0_yyzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 440);

            auto g_y_0_yyzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 441);

            auto g_y_0_yyzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 442);

            auto g_y_0_yyzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 443);

            auto g_y_0_yyzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 444);

            auto g_y_0_yyzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 445);

            auto g_y_0_yyzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 446);

            auto g_y_0_yyzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 447);

            auto g_y_0_yyzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 448);

            auto g_y_0_yyzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 449);

            auto g_y_0_yyzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 450);

            auto g_y_0_yyzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 451);

            auto g_y_0_yyzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 452);

            auto g_y_0_yyzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 453);

            auto g_y_0_yyzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 454);

            auto g_y_0_yyzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 455);

            auto g_y_0_yyzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 456);

            auto g_y_0_yyzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 457);

            auto g_y_0_yyzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 458);

            auto g_y_0_yyzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 459);

            auto g_y_0_yyzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 460);

            auto g_y_0_yyzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 461);

            auto g_y_0_yyzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 462);

            auto g_y_0_yyzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 463);

            auto g_y_0_yyzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 464);

            auto g_y_0_yyzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 465);

            auto g_y_0_yyzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 466);

            auto g_y_0_yyzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 467);

            auto g_y_0_yzzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 468);

            auto g_y_0_yzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 469);

            auto g_y_0_yzzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 470);

            auto g_y_0_yzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 471);

            auto g_y_0_yzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 472);

            auto g_y_0_yzzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 473);

            auto g_y_0_yzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 474);

            auto g_y_0_yzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 475);

            auto g_y_0_yzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 476);

            auto g_y_0_yzzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 477);

            auto g_y_0_yzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 478);

            auto g_y_0_yzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 479);

            auto g_y_0_yzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 480);

            auto g_y_0_yzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 481);

            auto g_y_0_yzzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 482);

            auto g_y_0_yzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 483);

            auto g_y_0_yzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 484);

            auto g_y_0_yzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 485);

            auto g_y_0_yzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 486);

            auto g_y_0_yzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 487);

            auto g_y_0_yzzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 488);

            auto g_y_0_yzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 489);

            auto g_y_0_yzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 490);

            auto g_y_0_yzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 491);

            auto g_y_0_yzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 492);

            auto g_y_0_yzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 493);

            auto g_y_0_yzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 494);

            auto g_y_0_yzzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 495);

            auto g_y_0_yzzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 496);

            auto g_y_0_yzzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 497);

            auto g_y_0_yzzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 498);

            auto g_y_0_yzzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 499);

            auto g_y_0_yzzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 500);

            auto g_y_0_yzzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 501);

            auto g_y_0_yzzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 502);

            auto g_y_0_yzzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 503);

            auto g_y_0_zzzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 504);

            auto g_y_0_zzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 505);

            auto g_y_0_zzzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 506);

            auto g_y_0_zzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 507);

            auto g_y_0_zzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 508);

            auto g_y_0_zzzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 509);

            auto g_y_0_zzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 510);

            auto g_y_0_zzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 511);

            auto g_y_0_zzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 512);

            auto g_y_0_zzzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 513);

            auto g_y_0_zzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 514);

            auto g_y_0_zzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 515);

            auto g_y_0_zzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 516);

            auto g_y_0_zzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 517);

            auto g_y_0_zzzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 518);

            auto g_y_0_zzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 519);

            auto g_y_0_zzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 520);

            auto g_y_0_zzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 521);

            auto g_y_0_zzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 522);

            auto g_y_0_zzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 523);

            auto g_y_0_zzzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 524);

            auto g_y_0_zzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 525);

            auto g_y_0_zzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 526);

            auto g_y_0_zzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 527);

            auto g_y_0_zzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 528);

            auto g_y_0_zzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 529);

            auto g_y_0_zzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 530);

            auto g_y_0_zzzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 531);

            auto g_y_0_zzzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 532);

            auto g_y_0_zzzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 533);

            auto g_y_0_zzzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 534);

            auto g_y_0_zzzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 535);

            auto g_y_0_zzzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 536);

            auto g_y_0_zzzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 537);

            auto g_y_0_zzzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 538);

            auto g_y_0_zzzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 540 * acomps * bcomps + 539);

            auto g_z_0_xxxx_xxxxxxx = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 0);

            auto g_z_0_xxxx_xxxxxxy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 1);

            auto g_z_0_xxxx_xxxxxxz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 2);

            auto g_z_0_xxxx_xxxxxyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 3);

            auto g_z_0_xxxx_xxxxxyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 4);

            auto g_z_0_xxxx_xxxxxzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 5);

            auto g_z_0_xxxx_xxxxyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 6);

            auto g_z_0_xxxx_xxxxyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 7);

            auto g_z_0_xxxx_xxxxyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 8);

            auto g_z_0_xxxx_xxxxzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 9);

            auto g_z_0_xxxx_xxxyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 10);

            auto g_z_0_xxxx_xxxyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 11);

            auto g_z_0_xxxx_xxxyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 12);

            auto g_z_0_xxxx_xxxyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 13);

            auto g_z_0_xxxx_xxxzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 14);

            auto g_z_0_xxxx_xxyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 15);

            auto g_z_0_xxxx_xxyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 16);

            auto g_z_0_xxxx_xxyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 17);

            auto g_z_0_xxxx_xxyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 18);

            auto g_z_0_xxxx_xxyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 19);

            auto g_z_0_xxxx_xxzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 20);

            auto g_z_0_xxxx_xyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 21);

            auto g_z_0_xxxx_xyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 22);

            auto g_z_0_xxxx_xyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 23);

            auto g_z_0_xxxx_xyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 24);

            auto g_z_0_xxxx_xyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 25);

            auto g_z_0_xxxx_xyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 26);

            auto g_z_0_xxxx_xzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 27);

            auto g_z_0_xxxx_yyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 28);

            auto g_z_0_xxxx_yyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 29);

            auto g_z_0_xxxx_yyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 30);

            auto g_z_0_xxxx_yyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 31);

            auto g_z_0_xxxx_yyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 32);

            auto g_z_0_xxxx_yyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 33);

            auto g_z_0_xxxx_yzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 34);

            auto g_z_0_xxxx_zzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 35);

            auto g_z_0_xxxy_xxxxxxx = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 36);

            auto g_z_0_xxxy_xxxxxxy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 37);

            auto g_z_0_xxxy_xxxxxxz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 38);

            auto g_z_0_xxxy_xxxxxyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 39);

            auto g_z_0_xxxy_xxxxxyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 40);

            auto g_z_0_xxxy_xxxxxzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 41);

            auto g_z_0_xxxy_xxxxyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 42);

            auto g_z_0_xxxy_xxxxyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 43);

            auto g_z_0_xxxy_xxxxyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 44);

            auto g_z_0_xxxy_xxxxzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 45);

            auto g_z_0_xxxy_xxxyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 46);

            auto g_z_0_xxxy_xxxyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 47);

            auto g_z_0_xxxy_xxxyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 48);

            auto g_z_0_xxxy_xxxyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 49);

            auto g_z_0_xxxy_xxxzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 50);

            auto g_z_0_xxxy_xxyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 51);

            auto g_z_0_xxxy_xxyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 52);

            auto g_z_0_xxxy_xxyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 53);

            auto g_z_0_xxxy_xxyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 54);

            auto g_z_0_xxxy_xxyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 55);

            auto g_z_0_xxxy_xxzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 56);

            auto g_z_0_xxxy_xyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 57);

            auto g_z_0_xxxy_xyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 58);

            auto g_z_0_xxxy_xyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 59);

            auto g_z_0_xxxy_xyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 60);

            auto g_z_0_xxxy_xyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 61);

            auto g_z_0_xxxy_xyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 62);

            auto g_z_0_xxxy_xzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 63);

            auto g_z_0_xxxy_yyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 64);

            auto g_z_0_xxxy_yyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 65);

            auto g_z_0_xxxy_yyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 66);

            auto g_z_0_xxxy_yyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 67);

            auto g_z_0_xxxy_yyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 68);

            auto g_z_0_xxxy_yyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 69);

            auto g_z_0_xxxy_yzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 70);

            auto g_z_0_xxxy_zzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 71);

            auto g_z_0_xxxz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 72);

            auto g_z_0_xxxz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 73);

            auto g_z_0_xxxz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 74);

            auto g_z_0_xxxz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 75);

            auto g_z_0_xxxz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 76);

            auto g_z_0_xxxz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 77);

            auto g_z_0_xxxz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 78);

            auto g_z_0_xxxz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 79);

            auto g_z_0_xxxz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 80);

            auto g_z_0_xxxz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 81);

            auto g_z_0_xxxz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 82);

            auto g_z_0_xxxz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 83);

            auto g_z_0_xxxz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 84);

            auto g_z_0_xxxz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 85);

            auto g_z_0_xxxz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 86);

            auto g_z_0_xxxz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 87);

            auto g_z_0_xxxz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 88);

            auto g_z_0_xxxz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 89);

            auto g_z_0_xxxz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 90);

            auto g_z_0_xxxz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 91);

            auto g_z_0_xxxz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 92);

            auto g_z_0_xxxz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 93);

            auto g_z_0_xxxz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 94);

            auto g_z_0_xxxz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 95);

            auto g_z_0_xxxz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 96);

            auto g_z_0_xxxz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 97);

            auto g_z_0_xxxz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 98);

            auto g_z_0_xxxz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 99);

            auto g_z_0_xxxz_yyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 100);

            auto g_z_0_xxxz_yyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 101);

            auto g_z_0_xxxz_yyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 102);

            auto g_z_0_xxxz_yyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 103);

            auto g_z_0_xxxz_yyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 104);

            auto g_z_0_xxxz_yyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 105);

            auto g_z_0_xxxz_yzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 106);

            auto g_z_0_xxxz_zzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 107);

            auto g_z_0_xxyy_xxxxxxx = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 108);

            auto g_z_0_xxyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 109);

            auto g_z_0_xxyy_xxxxxxz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 110);

            auto g_z_0_xxyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 111);

            auto g_z_0_xxyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 112);

            auto g_z_0_xxyy_xxxxxzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 113);

            auto g_z_0_xxyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 114);

            auto g_z_0_xxyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 115);

            auto g_z_0_xxyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 116);

            auto g_z_0_xxyy_xxxxzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 117);

            auto g_z_0_xxyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 118);

            auto g_z_0_xxyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 119);

            auto g_z_0_xxyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 120);

            auto g_z_0_xxyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 121);

            auto g_z_0_xxyy_xxxzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 122);

            auto g_z_0_xxyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 123);

            auto g_z_0_xxyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 124);

            auto g_z_0_xxyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 125);

            auto g_z_0_xxyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 126);

            auto g_z_0_xxyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 127);

            auto g_z_0_xxyy_xxzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 128);

            auto g_z_0_xxyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 129);

            auto g_z_0_xxyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 130);

            auto g_z_0_xxyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 131);

            auto g_z_0_xxyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 132);

            auto g_z_0_xxyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 133);

            auto g_z_0_xxyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 134);

            auto g_z_0_xxyy_xzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 135);

            auto g_z_0_xxyy_yyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 136);

            auto g_z_0_xxyy_yyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 137);

            auto g_z_0_xxyy_yyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 138);

            auto g_z_0_xxyy_yyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 139);

            auto g_z_0_xxyy_yyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 140);

            auto g_z_0_xxyy_yyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 141);

            auto g_z_0_xxyy_yzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 142);

            auto g_z_0_xxyy_zzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 143);

            auto g_z_0_xxyz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 144);

            auto g_z_0_xxyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 145);

            auto g_z_0_xxyz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 146);

            auto g_z_0_xxyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 147);

            auto g_z_0_xxyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 148);

            auto g_z_0_xxyz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 149);

            auto g_z_0_xxyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 150);

            auto g_z_0_xxyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 151);

            auto g_z_0_xxyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 152);

            auto g_z_0_xxyz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 153);

            auto g_z_0_xxyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 154);

            auto g_z_0_xxyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 155);

            auto g_z_0_xxyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 156);

            auto g_z_0_xxyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 157);

            auto g_z_0_xxyz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 158);

            auto g_z_0_xxyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 159);

            auto g_z_0_xxyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 160);

            auto g_z_0_xxyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 161);

            auto g_z_0_xxyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 162);

            auto g_z_0_xxyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 163);

            auto g_z_0_xxyz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 164);

            auto g_z_0_xxyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 165);

            auto g_z_0_xxyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 166);

            auto g_z_0_xxyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 167);

            auto g_z_0_xxyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 168);

            auto g_z_0_xxyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 169);

            auto g_z_0_xxyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 170);

            auto g_z_0_xxyz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 171);

            auto g_z_0_xxyz_yyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 172);

            auto g_z_0_xxyz_yyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 173);

            auto g_z_0_xxyz_yyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 174);

            auto g_z_0_xxyz_yyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 175);

            auto g_z_0_xxyz_yyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 176);

            auto g_z_0_xxyz_yyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 177);

            auto g_z_0_xxyz_yzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 178);

            auto g_z_0_xxyz_zzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 179);

            auto g_z_0_xxzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 180);

            auto g_z_0_xxzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 181);

            auto g_z_0_xxzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 182);

            auto g_z_0_xxzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 183);

            auto g_z_0_xxzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 184);

            auto g_z_0_xxzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 185);

            auto g_z_0_xxzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 186);

            auto g_z_0_xxzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 187);

            auto g_z_0_xxzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 188);

            auto g_z_0_xxzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 189);

            auto g_z_0_xxzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 190);

            auto g_z_0_xxzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 191);

            auto g_z_0_xxzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 192);

            auto g_z_0_xxzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 193);

            auto g_z_0_xxzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 194);

            auto g_z_0_xxzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 195);

            auto g_z_0_xxzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 196);

            auto g_z_0_xxzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 197);

            auto g_z_0_xxzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 198);

            auto g_z_0_xxzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 199);

            auto g_z_0_xxzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 200);

            auto g_z_0_xxzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 201);

            auto g_z_0_xxzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 202);

            auto g_z_0_xxzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 203);

            auto g_z_0_xxzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 204);

            auto g_z_0_xxzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 205);

            auto g_z_0_xxzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 206);

            auto g_z_0_xxzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 207);

            auto g_z_0_xxzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 208);

            auto g_z_0_xxzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 209);

            auto g_z_0_xxzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 210);

            auto g_z_0_xxzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 211);

            auto g_z_0_xxzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 212);

            auto g_z_0_xxzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 213);

            auto g_z_0_xxzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 214);

            auto g_z_0_xxzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 215);

            auto g_z_0_xyyy_xxxxxxx = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 216);

            auto g_z_0_xyyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 217);

            auto g_z_0_xyyy_xxxxxxz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 218);

            auto g_z_0_xyyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 219);

            auto g_z_0_xyyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 220);

            auto g_z_0_xyyy_xxxxxzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 221);

            auto g_z_0_xyyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 222);

            auto g_z_0_xyyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 223);

            auto g_z_0_xyyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 224);

            auto g_z_0_xyyy_xxxxzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 225);

            auto g_z_0_xyyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 226);

            auto g_z_0_xyyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 227);

            auto g_z_0_xyyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 228);

            auto g_z_0_xyyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 229);

            auto g_z_0_xyyy_xxxzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 230);

            auto g_z_0_xyyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 231);

            auto g_z_0_xyyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 232);

            auto g_z_0_xyyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 233);

            auto g_z_0_xyyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 234);

            auto g_z_0_xyyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 235);

            auto g_z_0_xyyy_xxzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 236);

            auto g_z_0_xyyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 237);

            auto g_z_0_xyyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 238);

            auto g_z_0_xyyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 239);

            auto g_z_0_xyyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 240);

            auto g_z_0_xyyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 241);

            auto g_z_0_xyyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 242);

            auto g_z_0_xyyy_xzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 243);

            auto g_z_0_xyyy_yyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 244);

            auto g_z_0_xyyy_yyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 245);

            auto g_z_0_xyyy_yyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 246);

            auto g_z_0_xyyy_yyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 247);

            auto g_z_0_xyyy_yyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 248);

            auto g_z_0_xyyy_yyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 249);

            auto g_z_0_xyyy_yzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 250);

            auto g_z_0_xyyy_zzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 251);

            auto g_z_0_xyyz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 252);

            auto g_z_0_xyyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 253);

            auto g_z_0_xyyz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 254);

            auto g_z_0_xyyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 255);

            auto g_z_0_xyyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 256);

            auto g_z_0_xyyz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 257);

            auto g_z_0_xyyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 258);

            auto g_z_0_xyyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 259);

            auto g_z_0_xyyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 260);

            auto g_z_0_xyyz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 261);

            auto g_z_0_xyyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 262);

            auto g_z_0_xyyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 263);

            auto g_z_0_xyyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 264);

            auto g_z_0_xyyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 265);

            auto g_z_0_xyyz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 266);

            auto g_z_0_xyyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 267);

            auto g_z_0_xyyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 268);

            auto g_z_0_xyyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 269);

            auto g_z_0_xyyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 270);

            auto g_z_0_xyyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 271);

            auto g_z_0_xyyz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 272);

            auto g_z_0_xyyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 273);

            auto g_z_0_xyyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 274);

            auto g_z_0_xyyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 275);

            auto g_z_0_xyyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 276);

            auto g_z_0_xyyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 277);

            auto g_z_0_xyyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 278);

            auto g_z_0_xyyz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 279);

            auto g_z_0_xyyz_yyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 280);

            auto g_z_0_xyyz_yyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 281);

            auto g_z_0_xyyz_yyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 282);

            auto g_z_0_xyyz_yyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 283);

            auto g_z_0_xyyz_yyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 284);

            auto g_z_0_xyyz_yyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 285);

            auto g_z_0_xyyz_yzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 286);

            auto g_z_0_xyyz_zzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 287);

            auto g_z_0_xyzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 288);

            auto g_z_0_xyzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 289);

            auto g_z_0_xyzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 290);

            auto g_z_0_xyzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 291);

            auto g_z_0_xyzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 292);

            auto g_z_0_xyzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 293);

            auto g_z_0_xyzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 294);

            auto g_z_0_xyzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 295);

            auto g_z_0_xyzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 296);

            auto g_z_0_xyzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 297);

            auto g_z_0_xyzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 298);

            auto g_z_0_xyzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 299);

            auto g_z_0_xyzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 300);

            auto g_z_0_xyzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 301);

            auto g_z_0_xyzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 302);

            auto g_z_0_xyzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 303);

            auto g_z_0_xyzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 304);

            auto g_z_0_xyzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 305);

            auto g_z_0_xyzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 306);

            auto g_z_0_xyzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 307);

            auto g_z_0_xyzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 308);

            auto g_z_0_xyzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 309);

            auto g_z_0_xyzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 310);

            auto g_z_0_xyzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 311);

            auto g_z_0_xyzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 312);

            auto g_z_0_xyzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 313);

            auto g_z_0_xyzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 314);

            auto g_z_0_xyzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 315);

            auto g_z_0_xyzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 316);

            auto g_z_0_xyzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 317);

            auto g_z_0_xyzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 318);

            auto g_z_0_xyzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 319);

            auto g_z_0_xyzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 320);

            auto g_z_0_xyzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 321);

            auto g_z_0_xyzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 322);

            auto g_z_0_xyzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 323);

            auto g_z_0_xzzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 324);

            auto g_z_0_xzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 325);

            auto g_z_0_xzzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 326);

            auto g_z_0_xzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 327);

            auto g_z_0_xzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 328);

            auto g_z_0_xzzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 329);

            auto g_z_0_xzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 330);

            auto g_z_0_xzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 331);

            auto g_z_0_xzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 332);

            auto g_z_0_xzzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 333);

            auto g_z_0_xzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 334);

            auto g_z_0_xzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 335);

            auto g_z_0_xzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 336);

            auto g_z_0_xzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 337);

            auto g_z_0_xzzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 338);

            auto g_z_0_xzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 339);

            auto g_z_0_xzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 340);

            auto g_z_0_xzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 341);

            auto g_z_0_xzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 342);

            auto g_z_0_xzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 343);

            auto g_z_0_xzzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 344);

            auto g_z_0_xzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 345);

            auto g_z_0_xzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 346);

            auto g_z_0_xzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 347);

            auto g_z_0_xzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 348);

            auto g_z_0_xzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 349);

            auto g_z_0_xzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 350);

            auto g_z_0_xzzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 351);

            auto g_z_0_xzzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 352);

            auto g_z_0_xzzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 353);

            auto g_z_0_xzzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 354);

            auto g_z_0_xzzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 355);

            auto g_z_0_xzzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 356);

            auto g_z_0_xzzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 357);

            auto g_z_0_xzzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 358);

            auto g_z_0_xzzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 359);

            auto g_z_0_yyyy_xxxxxxx = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 360);

            auto g_z_0_yyyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 361);

            auto g_z_0_yyyy_xxxxxxz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 362);

            auto g_z_0_yyyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 363);

            auto g_z_0_yyyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 364);

            auto g_z_0_yyyy_xxxxxzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 365);

            auto g_z_0_yyyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 366);

            auto g_z_0_yyyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 367);

            auto g_z_0_yyyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 368);

            auto g_z_0_yyyy_xxxxzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 369);

            auto g_z_0_yyyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 370);

            auto g_z_0_yyyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 371);

            auto g_z_0_yyyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 372);

            auto g_z_0_yyyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 373);

            auto g_z_0_yyyy_xxxzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 374);

            auto g_z_0_yyyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 375);

            auto g_z_0_yyyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 376);

            auto g_z_0_yyyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 377);

            auto g_z_0_yyyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 378);

            auto g_z_0_yyyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 379);

            auto g_z_0_yyyy_xxzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 380);

            auto g_z_0_yyyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 381);

            auto g_z_0_yyyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 382);

            auto g_z_0_yyyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 383);

            auto g_z_0_yyyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 384);

            auto g_z_0_yyyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 385);

            auto g_z_0_yyyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 386);

            auto g_z_0_yyyy_xzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 387);

            auto g_z_0_yyyy_yyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 388);

            auto g_z_0_yyyy_yyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 389);

            auto g_z_0_yyyy_yyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 390);

            auto g_z_0_yyyy_yyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 391);

            auto g_z_0_yyyy_yyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 392);

            auto g_z_0_yyyy_yyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 393);

            auto g_z_0_yyyy_yzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 394);

            auto g_z_0_yyyy_zzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 395);

            auto g_z_0_yyyz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 396);

            auto g_z_0_yyyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 397);

            auto g_z_0_yyyz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 398);

            auto g_z_0_yyyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 399);

            auto g_z_0_yyyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 400);

            auto g_z_0_yyyz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 401);

            auto g_z_0_yyyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 402);

            auto g_z_0_yyyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 403);

            auto g_z_0_yyyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 404);

            auto g_z_0_yyyz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 405);

            auto g_z_0_yyyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 406);

            auto g_z_0_yyyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 407);

            auto g_z_0_yyyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 408);

            auto g_z_0_yyyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 409);

            auto g_z_0_yyyz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 410);

            auto g_z_0_yyyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 411);

            auto g_z_0_yyyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 412);

            auto g_z_0_yyyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 413);

            auto g_z_0_yyyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 414);

            auto g_z_0_yyyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 415);

            auto g_z_0_yyyz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 416);

            auto g_z_0_yyyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 417);

            auto g_z_0_yyyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 418);

            auto g_z_0_yyyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 419);

            auto g_z_0_yyyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 420);

            auto g_z_0_yyyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 421);

            auto g_z_0_yyyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 422);

            auto g_z_0_yyyz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 423);

            auto g_z_0_yyyz_yyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 424);

            auto g_z_0_yyyz_yyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 425);

            auto g_z_0_yyyz_yyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 426);

            auto g_z_0_yyyz_yyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 427);

            auto g_z_0_yyyz_yyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 428);

            auto g_z_0_yyyz_yyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 429);

            auto g_z_0_yyyz_yzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 430);

            auto g_z_0_yyyz_zzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 431);

            auto g_z_0_yyzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 432);

            auto g_z_0_yyzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 433);

            auto g_z_0_yyzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 434);

            auto g_z_0_yyzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 435);

            auto g_z_0_yyzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 436);

            auto g_z_0_yyzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 437);

            auto g_z_0_yyzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 438);

            auto g_z_0_yyzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 439);

            auto g_z_0_yyzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 440);

            auto g_z_0_yyzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 441);

            auto g_z_0_yyzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 442);

            auto g_z_0_yyzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 443);

            auto g_z_0_yyzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 444);

            auto g_z_0_yyzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 445);

            auto g_z_0_yyzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 446);

            auto g_z_0_yyzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 447);

            auto g_z_0_yyzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 448);

            auto g_z_0_yyzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 449);

            auto g_z_0_yyzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 450);

            auto g_z_0_yyzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 451);

            auto g_z_0_yyzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 452);

            auto g_z_0_yyzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 453);

            auto g_z_0_yyzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 454);

            auto g_z_0_yyzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 455);

            auto g_z_0_yyzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 456);

            auto g_z_0_yyzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 457);

            auto g_z_0_yyzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 458);

            auto g_z_0_yyzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 459);

            auto g_z_0_yyzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 460);

            auto g_z_0_yyzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 461);

            auto g_z_0_yyzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 462);

            auto g_z_0_yyzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 463);

            auto g_z_0_yyzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 464);

            auto g_z_0_yyzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 465);

            auto g_z_0_yyzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 466);

            auto g_z_0_yyzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 467);

            auto g_z_0_yzzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 468);

            auto g_z_0_yzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 469);

            auto g_z_0_yzzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 470);

            auto g_z_0_yzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 471);

            auto g_z_0_yzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 472);

            auto g_z_0_yzzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 473);

            auto g_z_0_yzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 474);

            auto g_z_0_yzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 475);

            auto g_z_0_yzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 476);

            auto g_z_0_yzzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 477);

            auto g_z_0_yzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 478);

            auto g_z_0_yzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 479);

            auto g_z_0_yzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 480);

            auto g_z_0_yzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 481);

            auto g_z_0_yzzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 482);

            auto g_z_0_yzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 483);

            auto g_z_0_yzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 484);

            auto g_z_0_yzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 485);

            auto g_z_0_yzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 486);

            auto g_z_0_yzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 487);

            auto g_z_0_yzzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 488);

            auto g_z_0_yzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 489);

            auto g_z_0_yzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 490);

            auto g_z_0_yzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 491);

            auto g_z_0_yzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 492);

            auto g_z_0_yzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 493);

            auto g_z_0_yzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 494);

            auto g_z_0_yzzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 495);

            auto g_z_0_yzzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 496);

            auto g_z_0_yzzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 497);

            auto g_z_0_yzzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 498);

            auto g_z_0_yzzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 499);

            auto g_z_0_yzzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 500);

            auto g_z_0_yzzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 501);

            auto g_z_0_yzzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 502);

            auto g_z_0_yzzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 503);

            auto g_z_0_zzzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 504);

            auto g_z_0_zzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 505);

            auto g_z_0_zzzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 506);

            auto g_z_0_zzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 507);

            auto g_z_0_zzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 508);

            auto g_z_0_zzzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 509);

            auto g_z_0_zzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 510);

            auto g_z_0_zzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 511);

            auto g_z_0_zzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 512);

            auto g_z_0_zzzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 513);

            auto g_z_0_zzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 514);

            auto g_z_0_zzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 515);

            auto g_z_0_zzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 516);

            auto g_z_0_zzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 517);

            auto g_z_0_zzzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 518);

            auto g_z_0_zzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 519);

            auto g_z_0_zzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 520);

            auto g_z_0_zzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 521);

            auto g_z_0_zzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 522);

            auto g_z_0_zzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 523);

            auto g_z_0_zzzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 524);

            auto g_z_0_zzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 525);

            auto g_z_0_zzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 526);

            auto g_z_0_zzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 527);

            auto g_z_0_zzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 528);

            auto g_z_0_zzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 529);

            auto g_z_0_zzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 530);

            auto g_z_0_zzzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 531);

            auto g_z_0_zzzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 532);

            auto g_z_0_zzzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 533);

            auto g_z_0_zzzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 534);

            auto g_z_0_zzzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 535);

            auto g_z_0_zzzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 536);

            auto g_z_0_zzzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 537);

            auto g_z_0_zzzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 538);

            auto g_z_0_zzzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 1080 * acomps * bcomps + 539);

            /// set up bra offset for contr_buffer_xxhi

            const auto hi_geom_10_off = idx_geom_10_xxhi + (i * bcomps + j) * 588;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxx_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxx_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxx_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxxx_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxx_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxx_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxxx_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxx_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxx_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxxx_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxxx_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxxx_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxxx_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxxx_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxxx_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxxx_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxxx_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxxx_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxxx_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxxx_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxxx_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxxx_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxxx_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxxx_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxxx_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxxx_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxxx_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxxx_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 27);

            #pragma omp simd aligned(cd_x, g_x_0_xxxx_xxxxxx, g_x_0_xxxx_xxxxxxx, g_x_0_xxxx_xxxxxxy, g_x_0_xxxx_xxxxxxz, g_x_0_xxxx_xxxxxy, g_x_0_xxxx_xxxxxyy, g_x_0_xxxx_xxxxxyz, g_x_0_xxxx_xxxxxz, g_x_0_xxxx_xxxxxzz, g_x_0_xxxx_xxxxyy, g_x_0_xxxx_xxxxyyy, g_x_0_xxxx_xxxxyyz, g_x_0_xxxx_xxxxyz, g_x_0_xxxx_xxxxyzz, g_x_0_xxxx_xxxxzz, g_x_0_xxxx_xxxxzzz, g_x_0_xxxx_xxxyyy, g_x_0_xxxx_xxxyyyy, g_x_0_xxxx_xxxyyyz, g_x_0_xxxx_xxxyyz, g_x_0_xxxx_xxxyyzz, g_x_0_xxxx_xxxyzz, g_x_0_xxxx_xxxyzzz, g_x_0_xxxx_xxxzzz, g_x_0_xxxx_xxxzzzz, g_x_0_xxxx_xxyyyy, g_x_0_xxxx_xxyyyyy, g_x_0_xxxx_xxyyyyz, g_x_0_xxxx_xxyyyz, g_x_0_xxxx_xxyyyzz, g_x_0_xxxx_xxyyzz, g_x_0_xxxx_xxyyzzz, g_x_0_xxxx_xxyzzz, g_x_0_xxxx_xxyzzzz, g_x_0_xxxx_xxzzzz, g_x_0_xxxx_xxzzzzz, g_x_0_xxxx_xyyyyy, g_x_0_xxxx_xyyyyyy, g_x_0_xxxx_xyyyyyz, g_x_0_xxxx_xyyyyz, g_x_0_xxxx_xyyyyzz, g_x_0_xxxx_xyyyzz, g_x_0_xxxx_xyyyzzz, g_x_0_xxxx_xyyzzz, g_x_0_xxxx_xyyzzzz, g_x_0_xxxx_xyzzzz, g_x_0_xxxx_xyzzzzz, g_x_0_xxxx_xzzzzz, g_x_0_xxxx_xzzzzzz, g_x_0_xxxx_yyyyyy, g_x_0_xxxx_yyyyyz, g_x_0_xxxx_yyyyzz, g_x_0_xxxx_yyyzzz, g_x_0_xxxx_yyzzzz, g_x_0_xxxx_yzzzzz, g_x_0_xxxx_zzzzzz, g_x_0_xxxxx_xxxxxx, g_x_0_xxxxx_xxxxxy, g_x_0_xxxxx_xxxxxz, g_x_0_xxxxx_xxxxyy, g_x_0_xxxxx_xxxxyz, g_x_0_xxxxx_xxxxzz, g_x_0_xxxxx_xxxyyy, g_x_0_xxxxx_xxxyyz, g_x_0_xxxxx_xxxyzz, g_x_0_xxxxx_xxxzzz, g_x_0_xxxxx_xxyyyy, g_x_0_xxxxx_xxyyyz, g_x_0_xxxxx_xxyyzz, g_x_0_xxxxx_xxyzzz, g_x_0_xxxxx_xxzzzz, g_x_0_xxxxx_xyyyyy, g_x_0_xxxxx_xyyyyz, g_x_0_xxxxx_xyyyzz, g_x_0_xxxxx_xyyzzz, g_x_0_xxxxx_xyzzzz, g_x_0_xxxxx_xzzzzz, g_x_0_xxxxx_yyyyyy, g_x_0_xxxxx_yyyyyz, g_x_0_xxxxx_yyyyzz, g_x_0_xxxxx_yyyzzz, g_x_0_xxxxx_yyzzzz, g_x_0_xxxxx_yzzzzz, g_x_0_xxxxx_zzzzzz, g_xxxx_xxxxxx, g_xxxx_xxxxxy, g_xxxx_xxxxxz, g_xxxx_xxxxyy, g_xxxx_xxxxyz, g_xxxx_xxxxzz, g_xxxx_xxxyyy, g_xxxx_xxxyyz, g_xxxx_xxxyzz, g_xxxx_xxxzzz, g_xxxx_xxyyyy, g_xxxx_xxyyyz, g_xxxx_xxyyzz, g_xxxx_xxyzzz, g_xxxx_xxzzzz, g_xxxx_xyyyyy, g_xxxx_xyyyyz, g_xxxx_xyyyzz, g_xxxx_xyyzzz, g_xxxx_xyzzzz, g_xxxx_xzzzzz, g_xxxx_yyyyyy, g_xxxx_yyyyyz, g_xxxx_yyyyzz, g_xxxx_yyyzzz, g_xxxx_yyzzzz, g_xxxx_yzzzzz, g_xxxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxx_xxxxxx[k] = -g_xxxx_xxxxxx[k] - g_x_0_xxxx_xxxxxx[k] * cd_x[k] + g_x_0_xxxx_xxxxxxx[k];

                g_x_0_xxxxx_xxxxxy[k] = -g_xxxx_xxxxxy[k] - g_x_0_xxxx_xxxxxy[k] * cd_x[k] + g_x_0_xxxx_xxxxxxy[k];

                g_x_0_xxxxx_xxxxxz[k] = -g_xxxx_xxxxxz[k] - g_x_0_xxxx_xxxxxz[k] * cd_x[k] + g_x_0_xxxx_xxxxxxz[k];

                g_x_0_xxxxx_xxxxyy[k] = -g_xxxx_xxxxyy[k] - g_x_0_xxxx_xxxxyy[k] * cd_x[k] + g_x_0_xxxx_xxxxxyy[k];

                g_x_0_xxxxx_xxxxyz[k] = -g_xxxx_xxxxyz[k] - g_x_0_xxxx_xxxxyz[k] * cd_x[k] + g_x_0_xxxx_xxxxxyz[k];

                g_x_0_xxxxx_xxxxzz[k] = -g_xxxx_xxxxzz[k] - g_x_0_xxxx_xxxxzz[k] * cd_x[k] + g_x_0_xxxx_xxxxxzz[k];

                g_x_0_xxxxx_xxxyyy[k] = -g_xxxx_xxxyyy[k] - g_x_0_xxxx_xxxyyy[k] * cd_x[k] + g_x_0_xxxx_xxxxyyy[k];

                g_x_0_xxxxx_xxxyyz[k] = -g_xxxx_xxxyyz[k] - g_x_0_xxxx_xxxyyz[k] * cd_x[k] + g_x_0_xxxx_xxxxyyz[k];

                g_x_0_xxxxx_xxxyzz[k] = -g_xxxx_xxxyzz[k] - g_x_0_xxxx_xxxyzz[k] * cd_x[k] + g_x_0_xxxx_xxxxyzz[k];

                g_x_0_xxxxx_xxxzzz[k] = -g_xxxx_xxxzzz[k] - g_x_0_xxxx_xxxzzz[k] * cd_x[k] + g_x_0_xxxx_xxxxzzz[k];

                g_x_0_xxxxx_xxyyyy[k] = -g_xxxx_xxyyyy[k] - g_x_0_xxxx_xxyyyy[k] * cd_x[k] + g_x_0_xxxx_xxxyyyy[k];

                g_x_0_xxxxx_xxyyyz[k] = -g_xxxx_xxyyyz[k] - g_x_0_xxxx_xxyyyz[k] * cd_x[k] + g_x_0_xxxx_xxxyyyz[k];

                g_x_0_xxxxx_xxyyzz[k] = -g_xxxx_xxyyzz[k] - g_x_0_xxxx_xxyyzz[k] * cd_x[k] + g_x_0_xxxx_xxxyyzz[k];

                g_x_0_xxxxx_xxyzzz[k] = -g_xxxx_xxyzzz[k] - g_x_0_xxxx_xxyzzz[k] * cd_x[k] + g_x_0_xxxx_xxxyzzz[k];

                g_x_0_xxxxx_xxzzzz[k] = -g_xxxx_xxzzzz[k] - g_x_0_xxxx_xxzzzz[k] * cd_x[k] + g_x_0_xxxx_xxxzzzz[k];

                g_x_0_xxxxx_xyyyyy[k] = -g_xxxx_xyyyyy[k] - g_x_0_xxxx_xyyyyy[k] * cd_x[k] + g_x_0_xxxx_xxyyyyy[k];

                g_x_0_xxxxx_xyyyyz[k] = -g_xxxx_xyyyyz[k] - g_x_0_xxxx_xyyyyz[k] * cd_x[k] + g_x_0_xxxx_xxyyyyz[k];

                g_x_0_xxxxx_xyyyzz[k] = -g_xxxx_xyyyzz[k] - g_x_0_xxxx_xyyyzz[k] * cd_x[k] + g_x_0_xxxx_xxyyyzz[k];

                g_x_0_xxxxx_xyyzzz[k] = -g_xxxx_xyyzzz[k] - g_x_0_xxxx_xyyzzz[k] * cd_x[k] + g_x_0_xxxx_xxyyzzz[k];

                g_x_0_xxxxx_xyzzzz[k] = -g_xxxx_xyzzzz[k] - g_x_0_xxxx_xyzzzz[k] * cd_x[k] + g_x_0_xxxx_xxyzzzz[k];

                g_x_0_xxxxx_xzzzzz[k] = -g_xxxx_xzzzzz[k] - g_x_0_xxxx_xzzzzz[k] * cd_x[k] + g_x_0_xxxx_xxzzzzz[k];

                g_x_0_xxxxx_yyyyyy[k] = -g_xxxx_yyyyyy[k] - g_x_0_xxxx_yyyyyy[k] * cd_x[k] + g_x_0_xxxx_xyyyyyy[k];

                g_x_0_xxxxx_yyyyyz[k] = -g_xxxx_yyyyyz[k] - g_x_0_xxxx_yyyyyz[k] * cd_x[k] + g_x_0_xxxx_xyyyyyz[k];

                g_x_0_xxxxx_yyyyzz[k] = -g_xxxx_yyyyzz[k] - g_x_0_xxxx_yyyyzz[k] * cd_x[k] + g_x_0_xxxx_xyyyyzz[k];

                g_x_0_xxxxx_yyyzzz[k] = -g_xxxx_yyyzzz[k] - g_x_0_xxxx_yyyzzz[k] * cd_x[k] + g_x_0_xxxx_xyyyzzz[k];

                g_x_0_xxxxx_yyzzzz[k] = -g_xxxx_yyzzzz[k] - g_x_0_xxxx_yyzzzz[k] * cd_x[k] + g_x_0_xxxx_xyyzzzz[k];

                g_x_0_xxxxx_yzzzzz[k] = -g_xxxx_yzzzzz[k] - g_x_0_xxxx_yzzzzz[k] * cd_x[k] + g_x_0_xxxx_xyzzzzz[k];

                g_x_0_xxxxx_zzzzzz[k] = -g_xxxx_zzzzzz[k] - g_x_0_xxxx_zzzzzz[k] * cd_x[k] + g_x_0_xxxx_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxy_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxxy_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxxxy_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxxy_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxxy_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxxy_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxxy_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxxy_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxxy_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxxy_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxxy_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxxy_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxxy_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxxy_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxxxy_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxxy_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxxy_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxxxy_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxxy_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxxy_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxxxy_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxxy_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxxxy_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxxy_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxxy_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxxy_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxxy_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxxy_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 55);

            #pragma omp simd aligned(cd_y, g_x_0_xxxx_xxxxxx, g_x_0_xxxx_xxxxxxy, g_x_0_xxxx_xxxxxy, g_x_0_xxxx_xxxxxyy, g_x_0_xxxx_xxxxxyz, g_x_0_xxxx_xxxxxz, g_x_0_xxxx_xxxxyy, g_x_0_xxxx_xxxxyyy, g_x_0_xxxx_xxxxyyz, g_x_0_xxxx_xxxxyz, g_x_0_xxxx_xxxxyzz, g_x_0_xxxx_xxxxzz, g_x_0_xxxx_xxxyyy, g_x_0_xxxx_xxxyyyy, g_x_0_xxxx_xxxyyyz, g_x_0_xxxx_xxxyyz, g_x_0_xxxx_xxxyyzz, g_x_0_xxxx_xxxyzz, g_x_0_xxxx_xxxyzzz, g_x_0_xxxx_xxxzzz, g_x_0_xxxx_xxyyyy, g_x_0_xxxx_xxyyyyy, g_x_0_xxxx_xxyyyyz, g_x_0_xxxx_xxyyyz, g_x_0_xxxx_xxyyyzz, g_x_0_xxxx_xxyyzz, g_x_0_xxxx_xxyyzzz, g_x_0_xxxx_xxyzzz, g_x_0_xxxx_xxyzzzz, g_x_0_xxxx_xxzzzz, g_x_0_xxxx_xyyyyy, g_x_0_xxxx_xyyyyyy, g_x_0_xxxx_xyyyyyz, g_x_0_xxxx_xyyyyz, g_x_0_xxxx_xyyyyzz, g_x_0_xxxx_xyyyzz, g_x_0_xxxx_xyyyzzz, g_x_0_xxxx_xyyzzz, g_x_0_xxxx_xyyzzzz, g_x_0_xxxx_xyzzzz, g_x_0_xxxx_xyzzzzz, g_x_0_xxxx_xzzzzz, g_x_0_xxxx_yyyyyy, g_x_0_xxxx_yyyyyyy, g_x_0_xxxx_yyyyyyz, g_x_0_xxxx_yyyyyz, g_x_0_xxxx_yyyyyzz, g_x_0_xxxx_yyyyzz, g_x_0_xxxx_yyyyzzz, g_x_0_xxxx_yyyzzz, g_x_0_xxxx_yyyzzzz, g_x_0_xxxx_yyzzzz, g_x_0_xxxx_yyzzzzz, g_x_0_xxxx_yzzzzz, g_x_0_xxxx_yzzzzzz, g_x_0_xxxx_zzzzzz, g_x_0_xxxxy_xxxxxx, g_x_0_xxxxy_xxxxxy, g_x_0_xxxxy_xxxxxz, g_x_0_xxxxy_xxxxyy, g_x_0_xxxxy_xxxxyz, g_x_0_xxxxy_xxxxzz, g_x_0_xxxxy_xxxyyy, g_x_0_xxxxy_xxxyyz, g_x_0_xxxxy_xxxyzz, g_x_0_xxxxy_xxxzzz, g_x_0_xxxxy_xxyyyy, g_x_0_xxxxy_xxyyyz, g_x_0_xxxxy_xxyyzz, g_x_0_xxxxy_xxyzzz, g_x_0_xxxxy_xxzzzz, g_x_0_xxxxy_xyyyyy, g_x_0_xxxxy_xyyyyz, g_x_0_xxxxy_xyyyzz, g_x_0_xxxxy_xyyzzz, g_x_0_xxxxy_xyzzzz, g_x_0_xxxxy_xzzzzz, g_x_0_xxxxy_yyyyyy, g_x_0_xxxxy_yyyyyz, g_x_0_xxxxy_yyyyzz, g_x_0_xxxxy_yyyzzz, g_x_0_xxxxy_yyzzzz, g_x_0_xxxxy_yzzzzz, g_x_0_xxxxy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxy_xxxxxx[k] = -g_x_0_xxxx_xxxxxx[k] * cd_y[k] + g_x_0_xxxx_xxxxxxy[k];

                g_x_0_xxxxy_xxxxxy[k] = -g_x_0_xxxx_xxxxxy[k] * cd_y[k] + g_x_0_xxxx_xxxxxyy[k];

                g_x_0_xxxxy_xxxxxz[k] = -g_x_0_xxxx_xxxxxz[k] * cd_y[k] + g_x_0_xxxx_xxxxxyz[k];

                g_x_0_xxxxy_xxxxyy[k] = -g_x_0_xxxx_xxxxyy[k] * cd_y[k] + g_x_0_xxxx_xxxxyyy[k];

                g_x_0_xxxxy_xxxxyz[k] = -g_x_0_xxxx_xxxxyz[k] * cd_y[k] + g_x_0_xxxx_xxxxyyz[k];

                g_x_0_xxxxy_xxxxzz[k] = -g_x_0_xxxx_xxxxzz[k] * cd_y[k] + g_x_0_xxxx_xxxxyzz[k];

                g_x_0_xxxxy_xxxyyy[k] = -g_x_0_xxxx_xxxyyy[k] * cd_y[k] + g_x_0_xxxx_xxxyyyy[k];

                g_x_0_xxxxy_xxxyyz[k] = -g_x_0_xxxx_xxxyyz[k] * cd_y[k] + g_x_0_xxxx_xxxyyyz[k];

                g_x_0_xxxxy_xxxyzz[k] = -g_x_0_xxxx_xxxyzz[k] * cd_y[k] + g_x_0_xxxx_xxxyyzz[k];

                g_x_0_xxxxy_xxxzzz[k] = -g_x_0_xxxx_xxxzzz[k] * cd_y[k] + g_x_0_xxxx_xxxyzzz[k];

                g_x_0_xxxxy_xxyyyy[k] = -g_x_0_xxxx_xxyyyy[k] * cd_y[k] + g_x_0_xxxx_xxyyyyy[k];

                g_x_0_xxxxy_xxyyyz[k] = -g_x_0_xxxx_xxyyyz[k] * cd_y[k] + g_x_0_xxxx_xxyyyyz[k];

                g_x_0_xxxxy_xxyyzz[k] = -g_x_0_xxxx_xxyyzz[k] * cd_y[k] + g_x_0_xxxx_xxyyyzz[k];

                g_x_0_xxxxy_xxyzzz[k] = -g_x_0_xxxx_xxyzzz[k] * cd_y[k] + g_x_0_xxxx_xxyyzzz[k];

                g_x_0_xxxxy_xxzzzz[k] = -g_x_0_xxxx_xxzzzz[k] * cd_y[k] + g_x_0_xxxx_xxyzzzz[k];

                g_x_0_xxxxy_xyyyyy[k] = -g_x_0_xxxx_xyyyyy[k] * cd_y[k] + g_x_0_xxxx_xyyyyyy[k];

                g_x_0_xxxxy_xyyyyz[k] = -g_x_0_xxxx_xyyyyz[k] * cd_y[k] + g_x_0_xxxx_xyyyyyz[k];

                g_x_0_xxxxy_xyyyzz[k] = -g_x_0_xxxx_xyyyzz[k] * cd_y[k] + g_x_0_xxxx_xyyyyzz[k];

                g_x_0_xxxxy_xyyzzz[k] = -g_x_0_xxxx_xyyzzz[k] * cd_y[k] + g_x_0_xxxx_xyyyzzz[k];

                g_x_0_xxxxy_xyzzzz[k] = -g_x_0_xxxx_xyzzzz[k] * cd_y[k] + g_x_0_xxxx_xyyzzzz[k];

                g_x_0_xxxxy_xzzzzz[k] = -g_x_0_xxxx_xzzzzz[k] * cd_y[k] + g_x_0_xxxx_xyzzzzz[k];

                g_x_0_xxxxy_yyyyyy[k] = -g_x_0_xxxx_yyyyyy[k] * cd_y[k] + g_x_0_xxxx_yyyyyyy[k];

                g_x_0_xxxxy_yyyyyz[k] = -g_x_0_xxxx_yyyyyz[k] * cd_y[k] + g_x_0_xxxx_yyyyyyz[k];

                g_x_0_xxxxy_yyyyzz[k] = -g_x_0_xxxx_yyyyzz[k] * cd_y[k] + g_x_0_xxxx_yyyyyzz[k];

                g_x_0_xxxxy_yyyzzz[k] = -g_x_0_xxxx_yyyzzz[k] * cd_y[k] + g_x_0_xxxx_yyyyzzz[k];

                g_x_0_xxxxy_yyzzzz[k] = -g_x_0_xxxx_yyzzzz[k] * cd_y[k] + g_x_0_xxxx_yyyzzzz[k];

                g_x_0_xxxxy_yzzzzz[k] = -g_x_0_xxxx_yzzzzz[k] * cd_y[k] + g_x_0_xxxx_yyzzzzz[k];

                g_x_0_xxxxy_zzzzzz[k] = -g_x_0_xxxx_zzzzzz[k] * cd_y[k] + g_x_0_xxxx_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxxz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxxz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxxz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xxxxz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxxxz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxxxz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxxxz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxxxz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxxxz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxxxz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxxxz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxxxz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxxxz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxxxz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxxxz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxxxz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxxxz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxxxz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxxxz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxxxz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxxxz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxxxz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxxxz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxxxz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxxxz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxxxz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxxxz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_z, g_x_0_xxxx_xxxxxx, g_x_0_xxxx_xxxxxxz, g_x_0_xxxx_xxxxxy, g_x_0_xxxx_xxxxxyz, g_x_0_xxxx_xxxxxz, g_x_0_xxxx_xxxxxzz, g_x_0_xxxx_xxxxyy, g_x_0_xxxx_xxxxyyz, g_x_0_xxxx_xxxxyz, g_x_0_xxxx_xxxxyzz, g_x_0_xxxx_xxxxzz, g_x_0_xxxx_xxxxzzz, g_x_0_xxxx_xxxyyy, g_x_0_xxxx_xxxyyyz, g_x_0_xxxx_xxxyyz, g_x_0_xxxx_xxxyyzz, g_x_0_xxxx_xxxyzz, g_x_0_xxxx_xxxyzzz, g_x_0_xxxx_xxxzzz, g_x_0_xxxx_xxxzzzz, g_x_0_xxxx_xxyyyy, g_x_0_xxxx_xxyyyyz, g_x_0_xxxx_xxyyyz, g_x_0_xxxx_xxyyyzz, g_x_0_xxxx_xxyyzz, g_x_0_xxxx_xxyyzzz, g_x_0_xxxx_xxyzzz, g_x_0_xxxx_xxyzzzz, g_x_0_xxxx_xxzzzz, g_x_0_xxxx_xxzzzzz, g_x_0_xxxx_xyyyyy, g_x_0_xxxx_xyyyyyz, g_x_0_xxxx_xyyyyz, g_x_0_xxxx_xyyyyzz, g_x_0_xxxx_xyyyzz, g_x_0_xxxx_xyyyzzz, g_x_0_xxxx_xyyzzz, g_x_0_xxxx_xyyzzzz, g_x_0_xxxx_xyzzzz, g_x_0_xxxx_xyzzzzz, g_x_0_xxxx_xzzzzz, g_x_0_xxxx_xzzzzzz, g_x_0_xxxx_yyyyyy, g_x_0_xxxx_yyyyyyz, g_x_0_xxxx_yyyyyz, g_x_0_xxxx_yyyyyzz, g_x_0_xxxx_yyyyzz, g_x_0_xxxx_yyyyzzz, g_x_0_xxxx_yyyzzz, g_x_0_xxxx_yyyzzzz, g_x_0_xxxx_yyzzzz, g_x_0_xxxx_yyzzzzz, g_x_0_xxxx_yzzzzz, g_x_0_xxxx_yzzzzzz, g_x_0_xxxx_zzzzzz, g_x_0_xxxx_zzzzzzz, g_x_0_xxxxz_xxxxxx, g_x_0_xxxxz_xxxxxy, g_x_0_xxxxz_xxxxxz, g_x_0_xxxxz_xxxxyy, g_x_0_xxxxz_xxxxyz, g_x_0_xxxxz_xxxxzz, g_x_0_xxxxz_xxxyyy, g_x_0_xxxxz_xxxyyz, g_x_0_xxxxz_xxxyzz, g_x_0_xxxxz_xxxzzz, g_x_0_xxxxz_xxyyyy, g_x_0_xxxxz_xxyyyz, g_x_0_xxxxz_xxyyzz, g_x_0_xxxxz_xxyzzz, g_x_0_xxxxz_xxzzzz, g_x_0_xxxxz_xyyyyy, g_x_0_xxxxz_xyyyyz, g_x_0_xxxxz_xyyyzz, g_x_0_xxxxz_xyyzzz, g_x_0_xxxxz_xyzzzz, g_x_0_xxxxz_xzzzzz, g_x_0_xxxxz_yyyyyy, g_x_0_xxxxz_yyyyyz, g_x_0_xxxxz_yyyyzz, g_x_0_xxxxz_yyyzzz, g_x_0_xxxxz_yyzzzz, g_x_0_xxxxz_yzzzzz, g_x_0_xxxxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxz_xxxxxx[k] = -g_x_0_xxxx_xxxxxx[k] * cd_z[k] + g_x_0_xxxx_xxxxxxz[k];

                g_x_0_xxxxz_xxxxxy[k] = -g_x_0_xxxx_xxxxxy[k] * cd_z[k] + g_x_0_xxxx_xxxxxyz[k];

                g_x_0_xxxxz_xxxxxz[k] = -g_x_0_xxxx_xxxxxz[k] * cd_z[k] + g_x_0_xxxx_xxxxxzz[k];

                g_x_0_xxxxz_xxxxyy[k] = -g_x_0_xxxx_xxxxyy[k] * cd_z[k] + g_x_0_xxxx_xxxxyyz[k];

                g_x_0_xxxxz_xxxxyz[k] = -g_x_0_xxxx_xxxxyz[k] * cd_z[k] + g_x_0_xxxx_xxxxyzz[k];

                g_x_0_xxxxz_xxxxzz[k] = -g_x_0_xxxx_xxxxzz[k] * cd_z[k] + g_x_0_xxxx_xxxxzzz[k];

                g_x_0_xxxxz_xxxyyy[k] = -g_x_0_xxxx_xxxyyy[k] * cd_z[k] + g_x_0_xxxx_xxxyyyz[k];

                g_x_0_xxxxz_xxxyyz[k] = -g_x_0_xxxx_xxxyyz[k] * cd_z[k] + g_x_0_xxxx_xxxyyzz[k];

                g_x_0_xxxxz_xxxyzz[k] = -g_x_0_xxxx_xxxyzz[k] * cd_z[k] + g_x_0_xxxx_xxxyzzz[k];

                g_x_0_xxxxz_xxxzzz[k] = -g_x_0_xxxx_xxxzzz[k] * cd_z[k] + g_x_0_xxxx_xxxzzzz[k];

                g_x_0_xxxxz_xxyyyy[k] = -g_x_0_xxxx_xxyyyy[k] * cd_z[k] + g_x_0_xxxx_xxyyyyz[k];

                g_x_0_xxxxz_xxyyyz[k] = -g_x_0_xxxx_xxyyyz[k] * cd_z[k] + g_x_0_xxxx_xxyyyzz[k];

                g_x_0_xxxxz_xxyyzz[k] = -g_x_0_xxxx_xxyyzz[k] * cd_z[k] + g_x_0_xxxx_xxyyzzz[k];

                g_x_0_xxxxz_xxyzzz[k] = -g_x_0_xxxx_xxyzzz[k] * cd_z[k] + g_x_0_xxxx_xxyzzzz[k];

                g_x_0_xxxxz_xxzzzz[k] = -g_x_0_xxxx_xxzzzz[k] * cd_z[k] + g_x_0_xxxx_xxzzzzz[k];

                g_x_0_xxxxz_xyyyyy[k] = -g_x_0_xxxx_xyyyyy[k] * cd_z[k] + g_x_0_xxxx_xyyyyyz[k];

                g_x_0_xxxxz_xyyyyz[k] = -g_x_0_xxxx_xyyyyz[k] * cd_z[k] + g_x_0_xxxx_xyyyyzz[k];

                g_x_0_xxxxz_xyyyzz[k] = -g_x_0_xxxx_xyyyzz[k] * cd_z[k] + g_x_0_xxxx_xyyyzzz[k];

                g_x_0_xxxxz_xyyzzz[k] = -g_x_0_xxxx_xyyzzz[k] * cd_z[k] + g_x_0_xxxx_xyyzzzz[k];

                g_x_0_xxxxz_xyzzzz[k] = -g_x_0_xxxx_xyzzzz[k] * cd_z[k] + g_x_0_xxxx_xyzzzzz[k];

                g_x_0_xxxxz_xzzzzz[k] = -g_x_0_xxxx_xzzzzz[k] * cd_z[k] + g_x_0_xxxx_xzzzzzz[k];

                g_x_0_xxxxz_yyyyyy[k] = -g_x_0_xxxx_yyyyyy[k] * cd_z[k] + g_x_0_xxxx_yyyyyyz[k];

                g_x_0_xxxxz_yyyyyz[k] = -g_x_0_xxxx_yyyyyz[k] * cd_z[k] + g_x_0_xxxx_yyyyyzz[k];

                g_x_0_xxxxz_yyyyzz[k] = -g_x_0_xxxx_yyyyzz[k] * cd_z[k] + g_x_0_xxxx_yyyyzzz[k];

                g_x_0_xxxxz_yyyzzz[k] = -g_x_0_xxxx_yyyzzz[k] * cd_z[k] + g_x_0_xxxx_yyyzzzz[k];

                g_x_0_xxxxz_yyzzzz[k] = -g_x_0_xxxx_yyzzzz[k] * cd_z[k] + g_x_0_xxxx_yyzzzzz[k];

                g_x_0_xxxxz_yzzzzz[k] = -g_x_0_xxxx_yzzzzz[k] * cd_z[k] + g_x_0_xxxx_yzzzzzz[k];

                g_x_0_xxxxz_zzzzzz[k] = -g_x_0_xxxx_zzzzzz[k] * cd_z[k] + g_x_0_xxxx_zzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyy_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxxyy_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxxyy_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxxyy_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxxyy_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxxyy_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xxxyy_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxxyy_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxxyy_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxxyy_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxxyy_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxxyy_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxxyy_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxxyy_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxxyy_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxxyy_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xxxyy_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xxxyy_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xxxyy_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xxxyy_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xxxyy_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xxxyy_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xxxyy_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xxxyy_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xxxyy_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xxxyy_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xxxyy_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xxxyy_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 111);

            #pragma omp simd aligned(cd_y, g_x_0_xxxy_xxxxxx, g_x_0_xxxy_xxxxxxy, g_x_0_xxxy_xxxxxy, g_x_0_xxxy_xxxxxyy, g_x_0_xxxy_xxxxxyz, g_x_0_xxxy_xxxxxz, g_x_0_xxxy_xxxxyy, g_x_0_xxxy_xxxxyyy, g_x_0_xxxy_xxxxyyz, g_x_0_xxxy_xxxxyz, g_x_0_xxxy_xxxxyzz, g_x_0_xxxy_xxxxzz, g_x_0_xxxy_xxxyyy, g_x_0_xxxy_xxxyyyy, g_x_0_xxxy_xxxyyyz, g_x_0_xxxy_xxxyyz, g_x_0_xxxy_xxxyyzz, g_x_0_xxxy_xxxyzz, g_x_0_xxxy_xxxyzzz, g_x_0_xxxy_xxxzzz, g_x_0_xxxy_xxyyyy, g_x_0_xxxy_xxyyyyy, g_x_0_xxxy_xxyyyyz, g_x_0_xxxy_xxyyyz, g_x_0_xxxy_xxyyyzz, g_x_0_xxxy_xxyyzz, g_x_0_xxxy_xxyyzzz, g_x_0_xxxy_xxyzzz, g_x_0_xxxy_xxyzzzz, g_x_0_xxxy_xxzzzz, g_x_0_xxxy_xyyyyy, g_x_0_xxxy_xyyyyyy, g_x_0_xxxy_xyyyyyz, g_x_0_xxxy_xyyyyz, g_x_0_xxxy_xyyyyzz, g_x_0_xxxy_xyyyzz, g_x_0_xxxy_xyyyzzz, g_x_0_xxxy_xyyzzz, g_x_0_xxxy_xyyzzzz, g_x_0_xxxy_xyzzzz, g_x_0_xxxy_xyzzzzz, g_x_0_xxxy_xzzzzz, g_x_0_xxxy_yyyyyy, g_x_0_xxxy_yyyyyyy, g_x_0_xxxy_yyyyyyz, g_x_0_xxxy_yyyyyz, g_x_0_xxxy_yyyyyzz, g_x_0_xxxy_yyyyzz, g_x_0_xxxy_yyyyzzz, g_x_0_xxxy_yyyzzz, g_x_0_xxxy_yyyzzzz, g_x_0_xxxy_yyzzzz, g_x_0_xxxy_yyzzzzz, g_x_0_xxxy_yzzzzz, g_x_0_xxxy_yzzzzzz, g_x_0_xxxy_zzzzzz, g_x_0_xxxyy_xxxxxx, g_x_0_xxxyy_xxxxxy, g_x_0_xxxyy_xxxxxz, g_x_0_xxxyy_xxxxyy, g_x_0_xxxyy_xxxxyz, g_x_0_xxxyy_xxxxzz, g_x_0_xxxyy_xxxyyy, g_x_0_xxxyy_xxxyyz, g_x_0_xxxyy_xxxyzz, g_x_0_xxxyy_xxxzzz, g_x_0_xxxyy_xxyyyy, g_x_0_xxxyy_xxyyyz, g_x_0_xxxyy_xxyyzz, g_x_0_xxxyy_xxyzzz, g_x_0_xxxyy_xxzzzz, g_x_0_xxxyy_xyyyyy, g_x_0_xxxyy_xyyyyz, g_x_0_xxxyy_xyyyzz, g_x_0_xxxyy_xyyzzz, g_x_0_xxxyy_xyzzzz, g_x_0_xxxyy_xzzzzz, g_x_0_xxxyy_yyyyyy, g_x_0_xxxyy_yyyyyz, g_x_0_xxxyy_yyyyzz, g_x_0_xxxyy_yyyzzz, g_x_0_xxxyy_yyzzzz, g_x_0_xxxyy_yzzzzz, g_x_0_xxxyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyy_xxxxxx[k] = -g_x_0_xxxy_xxxxxx[k] * cd_y[k] + g_x_0_xxxy_xxxxxxy[k];

                g_x_0_xxxyy_xxxxxy[k] = -g_x_0_xxxy_xxxxxy[k] * cd_y[k] + g_x_0_xxxy_xxxxxyy[k];

                g_x_0_xxxyy_xxxxxz[k] = -g_x_0_xxxy_xxxxxz[k] * cd_y[k] + g_x_0_xxxy_xxxxxyz[k];

                g_x_0_xxxyy_xxxxyy[k] = -g_x_0_xxxy_xxxxyy[k] * cd_y[k] + g_x_0_xxxy_xxxxyyy[k];

                g_x_0_xxxyy_xxxxyz[k] = -g_x_0_xxxy_xxxxyz[k] * cd_y[k] + g_x_0_xxxy_xxxxyyz[k];

                g_x_0_xxxyy_xxxxzz[k] = -g_x_0_xxxy_xxxxzz[k] * cd_y[k] + g_x_0_xxxy_xxxxyzz[k];

                g_x_0_xxxyy_xxxyyy[k] = -g_x_0_xxxy_xxxyyy[k] * cd_y[k] + g_x_0_xxxy_xxxyyyy[k];

                g_x_0_xxxyy_xxxyyz[k] = -g_x_0_xxxy_xxxyyz[k] * cd_y[k] + g_x_0_xxxy_xxxyyyz[k];

                g_x_0_xxxyy_xxxyzz[k] = -g_x_0_xxxy_xxxyzz[k] * cd_y[k] + g_x_0_xxxy_xxxyyzz[k];

                g_x_0_xxxyy_xxxzzz[k] = -g_x_0_xxxy_xxxzzz[k] * cd_y[k] + g_x_0_xxxy_xxxyzzz[k];

                g_x_0_xxxyy_xxyyyy[k] = -g_x_0_xxxy_xxyyyy[k] * cd_y[k] + g_x_0_xxxy_xxyyyyy[k];

                g_x_0_xxxyy_xxyyyz[k] = -g_x_0_xxxy_xxyyyz[k] * cd_y[k] + g_x_0_xxxy_xxyyyyz[k];

                g_x_0_xxxyy_xxyyzz[k] = -g_x_0_xxxy_xxyyzz[k] * cd_y[k] + g_x_0_xxxy_xxyyyzz[k];

                g_x_0_xxxyy_xxyzzz[k] = -g_x_0_xxxy_xxyzzz[k] * cd_y[k] + g_x_0_xxxy_xxyyzzz[k];

                g_x_0_xxxyy_xxzzzz[k] = -g_x_0_xxxy_xxzzzz[k] * cd_y[k] + g_x_0_xxxy_xxyzzzz[k];

                g_x_0_xxxyy_xyyyyy[k] = -g_x_0_xxxy_xyyyyy[k] * cd_y[k] + g_x_0_xxxy_xyyyyyy[k];

                g_x_0_xxxyy_xyyyyz[k] = -g_x_0_xxxy_xyyyyz[k] * cd_y[k] + g_x_0_xxxy_xyyyyyz[k];

                g_x_0_xxxyy_xyyyzz[k] = -g_x_0_xxxy_xyyyzz[k] * cd_y[k] + g_x_0_xxxy_xyyyyzz[k];

                g_x_0_xxxyy_xyyzzz[k] = -g_x_0_xxxy_xyyzzz[k] * cd_y[k] + g_x_0_xxxy_xyyyzzz[k];

                g_x_0_xxxyy_xyzzzz[k] = -g_x_0_xxxy_xyzzzz[k] * cd_y[k] + g_x_0_xxxy_xyyzzzz[k];

                g_x_0_xxxyy_xzzzzz[k] = -g_x_0_xxxy_xzzzzz[k] * cd_y[k] + g_x_0_xxxy_xyzzzzz[k];

                g_x_0_xxxyy_yyyyyy[k] = -g_x_0_xxxy_yyyyyy[k] * cd_y[k] + g_x_0_xxxy_yyyyyyy[k];

                g_x_0_xxxyy_yyyyyz[k] = -g_x_0_xxxy_yyyyyz[k] * cd_y[k] + g_x_0_xxxy_yyyyyyz[k];

                g_x_0_xxxyy_yyyyzz[k] = -g_x_0_xxxy_yyyyzz[k] * cd_y[k] + g_x_0_xxxy_yyyyyzz[k];

                g_x_0_xxxyy_yyyzzz[k] = -g_x_0_xxxy_yyyzzz[k] * cd_y[k] + g_x_0_xxxy_yyyyzzz[k];

                g_x_0_xxxyy_yyzzzz[k] = -g_x_0_xxxy_yyzzzz[k] * cd_y[k] + g_x_0_xxxy_yyyzzzz[k];

                g_x_0_xxxyy_yzzzzz[k] = -g_x_0_xxxy_yzzzzz[k] * cd_y[k] + g_x_0_xxxy_yyzzzzz[k];

                g_x_0_xxxyy_zzzzzz[k] = -g_x_0_xxxy_zzzzzz[k] * cd_y[k] + g_x_0_xxxy_yzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xxxyz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xxxyz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xxxyz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xxxyz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xxxyz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xxxyz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xxxyz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xxxyz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xxxyz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xxxyz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xxxyz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xxxyz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xxxyz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xxxyz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xxxyz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xxxyz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xxxyz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xxxyz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xxxyz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xxxyz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xxxyz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xxxyz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xxxyz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xxxyz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xxxyz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xxxyz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xxxyz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 139);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyz_xxxxxx, g_x_0_xxxyz_xxxxxy, g_x_0_xxxyz_xxxxxz, g_x_0_xxxyz_xxxxyy, g_x_0_xxxyz_xxxxyz, g_x_0_xxxyz_xxxxzz, g_x_0_xxxyz_xxxyyy, g_x_0_xxxyz_xxxyyz, g_x_0_xxxyz_xxxyzz, g_x_0_xxxyz_xxxzzz, g_x_0_xxxyz_xxyyyy, g_x_0_xxxyz_xxyyyz, g_x_0_xxxyz_xxyyzz, g_x_0_xxxyz_xxyzzz, g_x_0_xxxyz_xxzzzz, g_x_0_xxxyz_xyyyyy, g_x_0_xxxyz_xyyyyz, g_x_0_xxxyz_xyyyzz, g_x_0_xxxyz_xyyzzz, g_x_0_xxxyz_xyzzzz, g_x_0_xxxyz_xzzzzz, g_x_0_xxxyz_yyyyyy, g_x_0_xxxyz_yyyyyz, g_x_0_xxxyz_yyyyzz, g_x_0_xxxyz_yyyzzz, g_x_0_xxxyz_yyzzzz, g_x_0_xxxyz_yzzzzz, g_x_0_xxxyz_zzzzzz, g_x_0_xxxz_xxxxxx, g_x_0_xxxz_xxxxxxy, g_x_0_xxxz_xxxxxy, g_x_0_xxxz_xxxxxyy, g_x_0_xxxz_xxxxxyz, g_x_0_xxxz_xxxxxz, g_x_0_xxxz_xxxxyy, g_x_0_xxxz_xxxxyyy, g_x_0_xxxz_xxxxyyz, g_x_0_xxxz_xxxxyz, g_x_0_xxxz_xxxxyzz, g_x_0_xxxz_xxxxzz, g_x_0_xxxz_xxxyyy, g_x_0_xxxz_xxxyyyy, g_x_0_xxxz_xxxyyyz, g_x_0_xxxz_xxxyyz, g_x_0_xxxz_xxxyyzz, g_x_0_xxxz_xxxyzz, g_x_0_xxxz_xxxyzzz, g_x_0_xxxz_xxxzzz, g_x_0_xxxz_xxyyyy, g_x_0_xxxz_xxyyyyy, g_x_0_xxxz_xxyyyyz, g_x_0_xxxz_xxyyyz, g_x_0_xxxz_xxyyyzz, g_x_0_xxxz_xxyyzz, g_x_0_xxxz_xxyyzzz, g_x_0_xxxz_xxyzzz, g_x_0_xxxz_xxyzzzz, g_x_0_xxxz_xxzzzz, g_x_0_xxxz_xyyyyy, g_x_0_xxxz_xyyyyyy, g_x_0_xxxz_xyyyyyz, g_x_0_xxxz_xyyyyz, g_x_0_xxxz_xyyyyzz, g_x_0_xxxz_xyyyzz, g_x_0_xxxz_xyyyzzz, g_x_0_xxxz_xyyzzz, g_x_0_xxxz_xyyzzzz, g_x_0_xxxz_xyzzzz, g_x_0_xxxz_xyzzzzz, g_x_0_xxxz_xzzzzz, g_x_0_xxxz_yyyyyy, g_x_0_xxxz_yyyyyyy, g_x_0_xxxz_yyyyyyz, g_x_0_xxxz_yyyyyz, g_x_0_xxxz_yyyyyzz, g_x_0_xxxz_yyyyzz, g_x_0_xxxz_yyyyzzz, g_x_0_xxxz_yyyzzz, g_x_0_xxxz_yyyzzzz, g_x_0_xxxz_yyzzzz, g_x_0_xxxz_yyzzzzz, g_x_0_xxxz_yzzzzz, g_x_0_xxxz_yzzzzzz, g_x_0_xxxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyz_xxxxxx[k] = -g_x_0_xxxz_xxxxxx[k] * cd_y[k] + g_x_0_xxxz_xxxxxxy[k];

                g_x_0_xxxyz_xxxxxy[k] = -g_x_0_xxxz_xxxxxy[k] * cd_y[k] + g_x_0_xxxz_xxxxxyy[k];

                g_x_0_xxxyz_xxxxxz[k] = -g_x_0_xxxz_xxxxxz[k] * cd_y[k] + g_x_0_xxxz_xxxxxyz[k];

                g_x_0_xxxyz_xxxxyy[k] = -g_x_0_xxxz_xxxxyy[k] * cd_y[k] + g_x_0_xxxz_xxxxyyy[k];

                g_x_0_xxxyz_xxxxyz[k] = -g_x_0_xxxz_xxxxyz[k] * cd_y[k] + g_x_0_xxxz_xxxxyyz[k];

                g_x_0_xxxyz_xxxxzz[k] = -g_x_0_xxxz_xxxxzz[k] * cd_y[k] + g_x_0_xxxz_xxxxyzz[k];

                g_x_0_xxxyz_xxxyyy[k] = -g_x_0_xxxz_xxxyyy[k] * cd_y[k] + g_x_0_xxxz_xxxyyyy[k];

                g_x_0_xxxyz_xxxyyz[k] = -g_x_0_xxxz_xxxyyz[k] * cd_y[k] + g_x_0_xxxz_xxxyyyz[k];

                g_x_0_xxxyz_xxxyzz[k] = -g_x_0_xxxz_xxxyzz[k] * cd_y[k] + g_x_0_xxxz_xxxyyzz[k];

                g_x_0_xxxyz_xxxzzz[k] = -g_x_0_xxxz_xxxzzz[k] * cd_y[k] + g_x_0_xxxz_xxxyzzz[k];

                g_x_0_xxxyz_xxyyyy[k] = -g_x_0_xxxz_xxyyyy[k] * cd_y[k] + g_x_0_xxxz_xxyyyyy[k];

                g_x_0_xxxyz_xxyyyz[k] = -g_x_0_xxxz_xxyyyz[k] * cd_y[k] + g_x_0_xxxz_xxyyyyz[k];

                g_x_0_xxxyz_xxyyzz[k] = -g_x_0_xxxz_xxyyzz[k] * cd_y[k] + g_x_0_xxxz_xxyyyzz[k];

                g_x_0_xxxyz_xxyzzz[k] = -g_x_0_xxxz_xxyzzz[k] * cd_y[k] + g_x_0_xxxz_xxyyzzz[k];

                g_x_0_xxxyz_xxzzzz[k] = -g_x_0_xxxz_xxzzzz[k] * cd_y[k] + g_x_0_xxxz_xxyzzzz[k];

                g_x_0_xxxyz_xyyyyy[k] = -g_x_0_xxxz_xyyyyy[k] * cd_y[k] + g_x_0_xxxz_xyyyyyy[k];

                g_x_0_xxxyz_xyyyyz[k] = -g_x_0_xxxz_xyyyyz[k] * cd_y[k] + g_x_0_xxxz_xyyyyyz[k];

                g_x_0_xxxyz_xyyyzz[k] = -g_x_0_xxxz_xyyyzz[k] * cd_y[k] + g_x_0_xxxz_xyyyyzz[k];

                g_x_0_xxxyz_xyyzzz[k] = -g_x_0_xxxz_xyyzzz[k] * cd_y[k] + g_x_0_xxxz_xyyyzzz[k];

                g_x_0_xxxyz_xyzzzz[k] = -g_x_0_xxxz_xyzzzz[k] * cd_y[k] + g_x_0_xxxz_xyyzzzz[k];

                g_x_0_xxxyz_xzzzzz[k] = -g_x_0_xxxz_xzzzzz[k] * cd_y[k] + g_x_0_xxxz_xyzzzzz[k];

                g_x_0_xxxyz_yyyyyy[k] = -g_x_0_xxxz_yyyyyy[k] * cd_y[k] + g_x_0_xxxz_yyyyyyy[k];

                g_x_0_xxxyz_yyyyyz[k] = -g_x_0_xxxz_yyyyyz[k] * cd_y[k] + g_x_0_xxxz_yyyyyyz[k];

                g_x_0_xxxyz_yyyyzz[k] = -g_x_0_xxxz_yyyyzz[k] * cd_y[k] + g_x_0_xxxz_yyyyyzz[k];

                g_x_0_xxxyz_yyyzzz[k] = -g_x_0_xxxz_yyyzzz[k] * cd_y[k] + g_x_0_xxxz_yyyyzzz[k];

                g_x_0_xxxyz_yyzzzz[k] = -g_x_0_xxxz_yyzzzz[k] * cd_y[k] + g_x_0_xxxz_yyyzzzz[k];

                g_x_0_xxxyz_yzzzzz[k] = -g_x_0_xxxz_yzzzzz[k] * cd_y[k] + g_x_0_xxxz_yyzzzzz[k];

                g_x_0_xxxyz_zzzzzz[k] = -g_x_0_xxxz_zzzzzz[k] * cd_y[k] + g_x_0_xxxz_yzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xxxzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xxxzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xxxzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xxxzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xxxzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xxxzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_xxxzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xxxzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xxxzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_xxxzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xxxzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xxxzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_xxxzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xxxzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xxxzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xxxzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_xxxzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xxxzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xxxzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_xxxzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xxxzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_xxxzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xxxzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xxxzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_xxxzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xxxzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xxxzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_z, g_x_0_xxxz_xxxxxx, g_x_0_xxxz_xxxxxxz, g_x_0_xxxz_xxxxxy, g_x_0_xxxz_xxxxxyz, g_x_0_xxxz_xxxxxz, g_x_0_xxxz_xxxxxzz, g_x_0_xxxz_xxxxyy, g_x_0_xxxz_xxxxyyz, g_x_0_xxxz_xxxxyz, g_x_0_xxxz_xxxxyzz, g_x_0_xxxz_xxxxzz, g_x_0_xxxz_xxxxzzz, g_x_0_xxxz_xxxyyy, g_x_0_xxxz_xxxyyyz, g_x_0_xxxz_xxxyyz, g_x_0_xxxz_xxxyyzz, g_x_0_xxxz_xxxyzz, g_x_0_xxxz_xxxyzzz, g_x_0_xxxz_xxxzzz, g_x_0_xxxz_xxxzzzz, g_x_0_xxxz_xxyyyy, g_x_0_xxxz_xxyyyyz, g_x_0_xxxz_xxyyyz, g_x_0_xxxz_xxyyyzz, g_x_0_xxxz_xxyyzz, g_x_0_xxxz_xxyyzzz, g_x_0_xxxz_xxyzzz, g_x_0_xxxz_xxyzzzz, g_x_0_xxxz_xxzzzz, g_x_0_xxxz_xxzzzzz, g_x_0_xxxz_xyyyyy, g_x_0_xxxz_xyyyyyz, g_x_0_xxxz_xyyyyz, g_x_0_xxxz_xyyyyzz, g_x_0_xxxz_xyyyzz, g_x_0_xxxz_xyyyzzz, g_x_0_xxxz_xyyzzz, g_x_0_xxxz_xyyzzzz, g_x_0_xxxz_xyzzzz, g_x_0_xxxz_xyzzzzz, g_x_0_xxxz_xzzzzz, g_x_0_xxxz_xzzzzzz, g_x_0_xxxz_yyyyyy, g_x_0_xxxz_yyyyyyz, g_x_0_xxxz_yyyyyz, g_x_0_xxxz_yyyyyzz, g_x_0_xxxz_yyyyzz, g_x_0_xxxz_yyyyzzz, g_x_0_xxxz_yyyzzz, g_x_0_xxxz_yyyzzzz, g_x_0_xxxz_yyzzzz, g_x_0_xxxz_yyzzzzz, g_x_0_xxxz_yzzzzz, g_x_0_xxxz_yzzzzzz, g_x_0_xxxz_zzzzzz, g_x_0_xxxz_zzzzzzz, g_x_0_xxxzz_xxxxxx, g_x_0_xxxzz_xxxxxy, g_x_0_xxxzz_xxxxxz, g_x_0_xxxzz_xxxxyy, g_x_0_xxxzz_xxxxyz, g_x_0_xxxzz_xxxxzz, g_x_0_xxxzz_xxxyyy, g_x_0_xxxzz_xxxyyz, g_x_0_xxxzz_xxxyzz, g_x_0_xxxzz_xxxzzz, g_x_0_xxxzz_xxyyyy, g_x_0_xxxzz_xxyyyz, g_x_0_xxxzz_xxyyzz, g_x_0_xxxzz_xxyzzz, g_x_0_xxxzz_xxzzzz, g_x_0_xxxzz_xyyyyy, g_x_0_xxxzz_xyyyyz, g_x_0_xxxzz_xyyyzz, g_x_0_xxxzz_xyyzzz, g_x_0_xxxzz_xyzzzz, g_x_0_xxxzz_xzzzzz, g_x_0_xxxzz_yyyyyy, g_x_0_xxxzz_yyyyyz, g_x_0_xxxzz_yyyyzz, g_x_0_xxxzz_yyyzzz, g_x_0_xxxzz_yyzzzz, g_x_0_xxxzz_yzzzzz, g_x_0_xxxzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzz_xxxxxx[k] = -g_x_0_xxxz_xxxxxx[k] * cd_z[k] + g_x_0_xxxz_xxxxxxz[k];

                g_x_0_xxxzz_xxxxxy[k] = -g_x_0_xxxz_xxxxxy[k] * cd_z[k] + g_x_0_xxxz_xxxxxyz[k];

                g_x_0_xxxzz_xxxxxz[k] = -g_x_0_xxxz_xxxxxz[k] * cd_z[k] + g_x_0_xxxz_xxxxxzz[k];

                g_x_0_xxxzz_xxxxyy[k] = -g_x_0_xxxz_xxxxyy[k] * cd_z[k] + g_x_0_xxxz_xxxxyyz[k];

                g_x_0_xxxzz_xxxxyz[k] = -g_x_0_xxxz_xxxxyz[k] * cd_z[k] + g_x_0_xxxz_xxxxyzz[k];

                g_x_0_xxxzz_xxxxzz[k] = -g_x_0_xxxz_xxxxzz[k] * cd_z[k] + g_x_0_xxxz_xxxxzzz[k];

                g_x_0_xxxzz_xxxyyy[k] = -g_x_0_xxxz_xxxyyy[k] * cd_z[k] + g_x_0_xxxz_xxxyyyz[k];

                g_x_0_xxxzz_xxxyyz[k] = -g_x_0_xxxz_xxxyyz[k] * cd_z[k] + g_x_0_xxxz_xxxyyzz[k];

                g_x_0_xxxzz_xxxyzz[k] = -g_x_0_xxxz_xxxyzz[k] * cd_z[k] + g_x_0_xxxz_xxxyzzz[k];

                g_x_0_xxxzz_xxxzzz[k] = -g_x_0_xxxz_xxxzzz[k] * cd_z[k] + g_x_0_xxxz_xxxzzzz[k];

                g_x_0_xxxzz_xxyyyy[k] = -g_x_0_xxxz_xxyyyy[k] * cd_z[k] + g_x_0_xxxz_xxyyyyz[k];

                g_x_0_xxxzz_xxyyyz[k] = -g_x_0_xxxz_xxyyyz[k] * cd_z[k] + g_x_0_xxxz_xxyyyzz[k];

                g_x_0_xxxzz_xxyyzz[k] = -g_x_0_xxxz_xxyyzz[k] * cd_z[k] + g_x_0_xxxz_xxyyzzz[k];

                g_x_0_xxxzz_xxyzzz[k] = -g_x_0_xxxz_xxyzzz[k] * cd_z[k] + g_x_0_xxxz_xxyzzzz[k];

                g_x_0_xxxzz_xxzzzz[k] = -g_x_0_xxxz_xxzzzz[k] * cd_z[k] + g_x_0_xxxz_xxzzzzz[k];

                g_x_0_xxxzz_xyyyyy[k] = -g_x_0_xxxz_xyyyyy[k] * cd_z[k] + g_x_0_xxxz_xyyyyyz[k];

                g_x_0_xxxzz_xyyyyz[k] = -g_x_0_xxxz_xyyyyz[k] * cd_z[k] + g_x_0_xxxz_xyyyyzz[k];

                g_x_0_xxxzz_xyyyzz[k] = -g_x_0_xxxz_xyyyzz[k] * cd_z[k] + g_x_0_xxxz_xyyyzzz[k];

                g_x_0_xxxzz_xyyzzz[k] = -g_x_0_xxxz_xyyzzz[k] * cd_z[k] + g_x_0_xxxz_xyyzzzz[k];

                g_x_0_xxxzz_xyzzzz[k] = -g_x_0_xxxz_xyzzzz[k] * cd_z[k] + g_x_0_xxxz_xyzzzzz[k];

                g_x_0_xxxzz_xzzzzz[k] = -g_x_0_xxxz_xzzzzz[k] * cd_z[k] + g_x_0_xxxz_xzzzzzz[k];

                g_x_0_xxxzz_yyyyyy[k] = -g_x_0_xxxz_yyyyyy[k] * cd_z[k] + g_x_0_xxxz_yyyyyyz[k];

                g_x_0_xxxzz_yyyyyz[k] = -g_x_0_xxxz_yyyyyz[k] * cd_z[k] + g_x_0_xxxz_yyyyyzz[k];

                g_x_0_xxxzz_yyyyzz[k] = -g_x_0_xxxz_yyyyzz[k] * cd_z[k] + g_x_0_xxxz_yyyyzzz[k];

                g_x_0_xxxzz_yyyzzz[k] = -g_x_0_xxxz_yyyzzz[k] * cd_z[k] + g_x_0_xxxz_yyyzzzz[k];

                g_x_0_xxxzz_yyzzzz[k] = -g_x_0_xxxz_yyzzzz[k] * cd_z[k] + g_x_0_xxxz_yyzzzzz[k];

                g_x_0_xxxzz_yzzzzz[k] = -g_x_0_xxxz_yzzzzz[k] * cd_z[k] + g_x_0_xxxz_yzzzzzz[k];

                g_x_0_xxxzz_zzzzzz[k] = -g_x_0_xxxz_zzzzzz[k] * cd_z[k] + g_x_0_xxxz_zzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_xxyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_xxyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_xxyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_xxyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_xxyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_xxyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_xxyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_xxyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_xxyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_xxyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_xxyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_xxyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_xxyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_xxyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_xxyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_xxyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_xxyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_xxyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_xxyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_xxyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_xxyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_xxyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_xxyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_xxyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_xxyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_xxyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_xxyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 195);

            #pragma omp simd aligned(cd_y, g_x_0_xxyy_xxxxxx, g_x_0_xxyy_xxxxxxy, g_x_0_xxyy_xxxxxy, g_x_0_xxyy_xxxxxyy, g_x_0_xxyy_xxxxxyz, g_x_0_xxyy_xxxxxz, g_x_0_xxyy_xxxxyy, g_x_0_xxyy_xxxxyyy, g_x_0_xxyy_xxxxyyz, g_x_0_xxyy_xxxxyz, g_x_0_xxyy_xxxxyzz, g_x_0_xxyy_xxxxzz, g_x_0_xxyy_xxxyyy, g_x_0_xxyy_xxxyyyy, g_x_0_xxyy_xxxyyyz, g_x_0_xxyy_xxxyyz, g_x_0_xxyy_xxxyyzz, g_x_0_xxyy_xxxyzz, g_x_0_xxyy_xxxyzzz, g_x_0_xxyy_xxxzzz, g_x_0_xxyy_xxyyyy, g_x_0_xxyy_xxyyyyy, g_x_0_xxyy_xxyyyyz, g_x_0_xxyy_xxyyyz, g_x_0_xxyy_xxyyyzz, g_x_0_xxyy_xxyyzz, g_x_0_xxyy_xxyyzzz, g_x_0_xxyy_xxyzzz, g_x_0_xxyy_xxyzzzz, g_x_0_xxyy_xxzzzz, g_x_0_xxyy_xyyyyy, g_x_0_xxyy_xyyyyyy, g_x_0_xxyy_xyyyyyz, g_x_0_xxyy_xyyyyz, g_x_0_xxyy_xyyyyzz, g_x_0_xxyy_xyyyzz, g_x_0_xxyy_xyyyzzz, g_x_0_xxyy_xyyzzz, g_x_0_xxyy_xyyzzzz, g_x_0_xxyy_xyzzzz, g_x_0_xxyy_xyzzzzz, g_x_0_xxyy_xzzzzz, g_x_0_xxyy_yyyyyy, g_x_0_xxyy_yyyyyyy, g_x_0_xxyy_yyyyyyz, g_x_0_xxyy_yyyyyz, g_x_0_xxyy_yyyyyzz, g_x_0_xxyy_yyyyzz, g_x_0_xxyy_yyyyzzz, g_x_0_xxyy_yyyzzz, g_x_0_xxyy_yyyzzzz, g_x_0_xxyy_yyzzzz, g_x_0_xxyy_yyzzzzz, g_x_0_xxyy_yzzzzz, g_x_0_xxyy_yzzzzzz, g_x_0_xxyy_zzzzzz, g_x_0_xxyyy_xxxxxx, g_x_0_xxyyy_xxxxxy, g_x_0_xxyyy_xxxxxz, g_x_0_xxyyy_xxxxyy, g_x_0_xxyyy_xxxxyz, g_x_0_xxyyy_xxxxzz, g_x_0_xxyyy_xxxyyy, g_x_0_xxyyy_xxxyyz, g_x_0_xxyyy_xxxyzz, g_x_0_xxyyy_xxxzzz, g_x_0_xxyyy_xxyyyy, g_x_0_xxyyy_xxyyyz, g_x_0_xxyyy_xxyyzz, g_x_0_xxyyy_xxyzzz, g_x_0_xxyyy_xxzzzz, g_x_0_xxyyy_xyyyyy, g_x_0_xxyyy_xyyyyz, g_x_0_xxyyy_xyyyzz, g_x_0_xxyyy_xyyzzz, g_x_0_xxyyy_xyzzzz, g_x_0_xxyyy_xzzzzz, g_x_0_xxyyy_yyyyyy, g_x_0_xxyyy_yyyyyz, g_x_0_xxyyy_yyyyzz, g_x_0_xxyyy_yyyzzz, g_x_0_xxyyy_yyzzzz, g_x_0_xxyyy_yzzzzz, g_x_0_xxyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyy_xxxxxx[k] = -g_x_0_xxyy_xxxxxx[k] * cd_y[k] + g_x_0_xxyy_xxxxxxy[k];

                g_x_0_xxyyy_xxxxxy[k] = -g_x_0_xxyy_xxxxxy[k] * cd_y[k] + g_x_0_xxyy_xxxxxyy[k];

                g_x_0_xxyyy_xxxxxz[k] = -g_x_0_xxyy_xxxxxz[k] * cd_y[k] + g_x_0_xxyy_xxxxxyz[k];

                g_x_0_xxyyy_xxxxyy[k] = -g_x_0_xxyy_xxxxyy[k] * cd_y[k] + g_x_0_xxyy_xxxxyyy[k];

                g_x_0_xxyyy_xxxxyz[k] = -g_x_0_xxyy_xxxxyz[k] * cd_y[k] + g_x_0_xxyy_xxxxyyz[k];

                g_x_0_xxyyy_xxxxzz[k] = -g_x_0_xxyy_xxxxzz[k] * cd_y[k] + g_x_0_xxyy_xxxxyzz[k];

                g_x_0_xxyyy_xxxyyy[k] = -g_x_0_xxyy_xxxyyy[k] * cd_y[k] + g_x_0_xxyy_xxxyyyy[k];

                g_x_0_xxyyy_xxxyyz[k] = -g_x_0_xxyy_xxxyyz[k] * cd_y[k] + g_x_0_xxyy_xxxyyyz[k];

                g_x_0_xxyyy_xxxyzz[k] = -g_x_0_xxyy_xxxyzz[k] * cd_y[k] + g_x_0_xxyy_xxxyyzz[k];

                g_x_0_xxyyy_xxxzzz[k] = -g_x_0_xxyy_xxxzzz[k] * cd_y[k] + g_x_0_xxyy_xxxyzzz[k];

                g_x_0_xxyyy_xxyyyy[k] = -g_x_0_xxyy_xxyyyy[k] * cd_y[k] + g_x_0_xxyy_xxyyyyy[k];

                g_x_0_xxyyy_xxyyyz[k] = -g_x_0_xxyy_xxyyyz[k] * cd_y[k] + g_x_0_xxyy_xxyyyyz[k];

                g_x_0_xxyyy_xxyyzz[k] = -g_x_0_xxyy_xxyyzz[k] * cd_y[k] + g_x_0_xxyy_xxyyyzz[k];

                g_x_0_xxyyy_xxyzzz[k] = -g_x_0_xxyy_xxyzzz[k] * cd_y[k] + g_x_0_xxyy_xxyyzzz[k];

                g_x_0_xxyyy_xxzzzz[k] = -g_x_0_xxyy_xxzzzz[k] * cd_y[k] + g_x_0_xxyy_xxyzzzz[k];

                g_x_0_xxyyy_xyyyyy[k] = -g_x_0_xxyy_xyyyyy[k] * cd_y[k] + g_x_0_xxyy_xyyyyyy[k];

                g_x_0_xxyyy_xyyyyz[k] = -g_x_0_xxyy_xyyyyz[k] * cd_y[k] + g_x_0_xxyy_xyyyyyz[k];

                g_x_0_xxyyy_xyyyzz[k] = -g_x_0_xxyy_xyyyzz[k] * cd_y[k] + g_x_0_xxyy_xyyyyzz[k];

                g_x_0_xxyyy_xyyzzz[k] = -g_x_0_xxyy_xyyzzz[k] * cd_y[k] + g_x_0_xxyy_xyyyzzz[k];

                g_x_0_xxyyy_xyzzzz[k] = -g_x_0_xxyy_xyzzzz[k] * cd_y[k] + g_x_0_xxyy_xyyzzzz[k];

                g_x_0_xxyyy_xzzzzz[k] = -g_x_0_xxyy_xzzzzz[k] * cd_y[k] + g_x_0_xxyy_xyzzzzz[k];

                g_x_0_xxyyy_yyyyyy[k] = -g_x_0_xxyy_yyyyyy[k] * cd_y[k] + g_x_0_xxyy_yyyyyyy[k];

                g_x_0_xxyyy_yyyyyz[k] = -g_x_0_xxyy_yyyyyz[k] * cd_y[k] + g_x_0_xxyy_yyyyyyz[k];

                g_x_0_xxyyy_yyyyzz[k] = -g_x_0_xxyy_yyyyzz[k] * cd_y[k] + g_x_0_xxyy_yyyyyzz[k];

                g_x_0_xxyyy_yyyzzz[k] = -g_x_0_xxyy_yyyzzz[k] * cd_y[k] + g_x_0_xxyy_yyyyzzz[k];

                g_x_0_xxyyy_yyzzzz[k] = -g_x_0_xxyy_yyzzzz[k] * cd_y[k] + g_x_0_xxyy_yyyzzzz[k];

                g_x_0_xxyyy_yzzzzz[k] = -g_x_0_xxyy_yzzzzz[k] * cd_y[k] + g_x_0_xxyy_yyzzzzz[k];

                g_x_0_xxyyy_zzzzzz[k] = -g_x_0_xxyy_zzzzzz[k] * cd_y[k] + g_x_0_xxyy_yzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_xxyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_xxyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_xxyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_xxyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_xxyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_xxyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_xxyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_xxyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_xxyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_xxyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_xxyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_xxyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_xxyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 209);

            auto g_x_0_xxyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_xxyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_xxyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_xxyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_xxyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_xxyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_xxyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_xxyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_xxyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_xxyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_xxyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_xxyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_xxyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_xxyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 223);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyz_xxxxxx, g_x_0_xxyyz_xxxxxy, g_x_0_xxyyz_xxxxxz, g_x_0_xxyyz_xxxxyy, g_x_0_xxyyz_xxxxyz, g_x_0_xxyyz_xxxxzz, g_x_0_xxyyz_xxxyyy, g_x_0_xxyyz_xxxyyz, g_x_0_xxyyz_xxxyzz, g_x_0_xxyyz_xxxzzz, g_x_0_xxyyz_xxyyyy, g_x_0_xxyyz_xxyyyz, g_x_0_xxyyz_xxyyzz, g_x_0_xxyyz_xxyzzz, g_x_0_xxyyz_xxzzzz, g_x_0_xxyyz_xyyyyy, g_x_0_xxyyz_xyyyyz, g_x_0_xxyyz_xyyyzz, g_x_0_xxyyz_xyyzzz, g_x_0_xxyyz_xyzzzz, g_x_0_xxyyz_xzzzzz, g_x_0_xxyyz_yyyyyy, g_x_0_xxyyz_yyyyyz, g_x_0_xxyyz_yyyyzz, g_x_0_xxyyz_yyyzzz, g_x_0_xxyyz_yyzzzz, g_x_0_xxyyz_yzzzzz, g_x_0_xxyyz_zzzzzz, g_x_0_xxyz_xxxxxx, g_x_0_xxyz_xxxxxxy, g_x_0_xxyz_xxxxxy, g_x_0_xxyz_xxxxxyy, g_x_0_xxyz_xxxxxyz, g_x_0_xxyz_xxxxxz, g_x_0_xxyz_xxxxyy, g_x_0_xxyz_xxxxyyy, g_x_0_xxyz_xxxxyyz, g_x_0_xxyz_xxxxyz, g_x_0_xxyz_xxxxyzz, g_x_0_xxyz_xxxxzz, g_x_0_xxyz_xxxyyy, g_x_0_xxyz_xxxyyyy, g_x_0_xxyz_xxxyyyz, g_x_0_xxyz_xxxyyz, g_x_0_xxyz_xxxyyzz, g_x_0_xxyz_xxxyzz, g_x_0_xxyz_xxxyzzz, g_x_0_xxyz_xxxzzz, g_x_0_xxyz_xxyyyy, g_x_0_xxyz_xxyyyyy, g_x_0_xxyz_xxyyyyz, g_x_0_xxyz_xxyyyz, g_x_0_xxyz_xxyyyzz, g_x_0_xxyz_xxyyzz, g_x_0_xxyz_xxyyzzz, g_x_0_xxyz_xxyzzz, g_x_0_xxyz_xxyzzzz, g_x_0_xxyz_xxzzzz, g_x_0_xxyz_xyyyyy, g_x_0_xxyz_xyyyyyy, g_x_0_xxyz_xyyyyyz, g_x_0_xxyz_xyyyyz, g_x_0_xxyz_xyyyyzz, g_x_0_xxyz_xyyyzz, g_x_0_xxyz_xyyyzzz, g_x_0_xxyz_xyyzzz, g_x_0_xxyz_xyyzzzz, g_x_0_xxyz_xyzzzz, g_x_0_xxyz_xyzzzzz, g_x_0_xxyz_xzzzzz, g_x_0_xxyz_yyyyyy, g_x_0_xxyz_yyyyyyy, g_x_0_xxyz_yyyyyyz, g_x_0_xxyz_yyyyyz, g_x_0_xxyz_yyyyyzz, g_x_0_xxyz_yyyyzz, g_x_0_xxyz_yyyyzzz, g_x_0_xxyz_yyyzzz, g_x_0_xxyz_yyyzzzz, g_x_0_xxyz_yyzzzz, g_x_0_xxyz_yyzzzzz, g_x_0_xxyz_yzzzzz, g_x_0_xxyz_yzzzzzz, g_x_0_xxyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyz_xxxxxx[k] = -g_x_0_xxyz_xxxxxx[k] * cd_y[k] + g_x_0_xxyz_xxxxxxy[k];

                g_x_0_xxyyz_xxxxxy[k] = -g_x_0_xxyz_xxxxxy[k] * cd_y[k] + g_x_0_xxyz_xxxxxyy[k];

                g_x_0_xxyyz_xxxxxz[k] = -g_x_0_xxyz_xxxxxz[k] * cd_y[k] + g_x_0_xxyz_xxxxxyz[k];

                g_x_0_xxyyz_xxxxyy[k] = -g_x_0_xxyz_xxxxyy[k] * cd_y[k] + g_x_0_xxyz_xxxxyyy[k];

                g_x_0_xxyyz_xxxxyz[k] = -g_x_0_xxyz_xxxxyz[k] * cd_y[k] + g_x_0_xxyz_xxxxyyz[k];

                g_x_0_xxyyz_xxxxzz[k] = -g_x_0_xxyz_xxxxzz[k] * cd_y[k] + g_x_0_xxyz_xxxxyzz[k];

                g_x_0_xxyyz_xxxyyy[k] = -g_x_0_xxyz_xxxyyy[k] * cd_y[k] + g_x_0_xxyz_xxxyyyy[k];

                g_x_0_xxyyz_xxxyyz[k] = -g_x_0_xxyz_xxxyyz[k] * cd_y[k] + g_x_0_xxyz_xxxyyyz[k];

                g_x_0_xxyyz_xxxyzz[k] = -g_x_0_xxyz_xxxyzz[k] * cd_y[k] + g_x_0_xxyz_xxxyyzz[k];

                g_x_0_xxyyz_xxxzzz[k] = -g_x_0_xxyz_xxxzzz[k] * cd_y[k] + g_x_0_xxyz_xxxyzzz[k];

                g_x_0_xxyyz_xxyyyy[k] = -g_x_0_xxyz_xxyyyy[k] * cd_y[k] + g_x_0_xxyz_xxyyyyy[k];

                g_x_0_xxyyz_xxyyyz[k] = -g_x_0_xxyz_xxyyyz[k] * cd_y[k] + g_x_0_xxyz_xxyyyyz[k];

                g_x_0_xxyyz_xxyyzz[k] = -g_x_0_xxyz_xxyyzz[k] * cd_y[k] + g_x_0_xxyz_xxyyyzz[k];

                g_x_0_xxyyz_xxyzzz[k] = -g_x_0_xxyz_xxyzzz[k] * cd_y[k] + g_x_0_xxyz_xxyyzzz[k];

                g_x_0_xxyyz_xxzzzz[k] = -g_x_0_xxyz_xxzzzz[k] * cd_y[k] + g_x_0_xxyz_xxyzzzz[k];

                g_x_0_xxyyz_xyyyyy[k] = -g_x_0_xxyz_xyyyyy[k] * cd_y[k] + g_x_0_xxyz_xyyyyyy[k];

                g_x_0_xxyyz_xyyyyz[k] = -g_x_0_xxyz_xyyyyz[k] * cd_y[k] + g_x_0_xxyz_xyyyyyz[k];

                g_x_0_xxyyz_xyyyzz[k] = -g_x_0_xxyz_xyyyzz[k] * cd_y[k] + g_x_0_xxyz_xyyyyzz[k];

                g_x_0_xxyyz_xyyzzz[k] = -g_x_0_xxyz_xyyzzz[k] * cd_y[k] + g_x_0_xxyz_xyyyzzz[k];

                g_x_0_xxyyz_xyzzzz[k] = -g_x_0_xxyz_xyzzzz[k] * cd_y[k] + g_x_0_xxyz_xyyzzzz[k];

                g_x_0_xxyyz_xzzzzz[k] = -g_x_0_xxyz_xzzzzz[k] * cd_y[k] + g_x_0_xxyz_xyzzzzz[k];

                g_x_0_xxyyz_yyyyyy[k] = -g_x_0_xxyz_yyyyyy[k] * cd_y[k] + g_x_0_xxyz_yyyyyyy[k];

                g_x_0_xxyyz_yyyyyz[k] = -g_x_0_xxyz_yyyyyz[k] * cd_y[k] + g_x_0_xxyz_yyyyyyz[k];

                g_x_0_xxyyz_yyyyzz[k] = -g_x_0_xxyz_yyyyzz[k] * cd_y[k] + g_x_0_xxyz_yyyyyzz[k];

                g_x_0_xxyyz_yyyzzz[k] = -g_x_0_xxyz_yyyzzz[k] * cd_y[k] + g_x_0_xxyz_yyyyzzz[k];

                g_x_0_xxyyz_yyzzzz[k] = -g_x_0_xxyz_yyzzzz[k] * cd_y[k] + g_x_0_xxyz_yyyzzzz[k];

                g_x_0_xxyyz_yzzzzz[k] = -g_x_0_xxyz_yzzzzz[k] * cd_y[k] + g_x_0_xxyz_yyzzzzz[k];

                g_x_0_xxyyz_zzzzzz[k] = -g_x_0_xxyz_zzzzzz[k] * cd_y[k] + g_x_0_xxyz_yzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 224);

            auto g_x_0_xxyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_xxyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_xxyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_xxyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_xxyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 229);

            auto g_x_0_xxyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 230);

            auto g_x_0_xxyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_xxyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_xxyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 233);

            auto g_x_0_xxyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_xxyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_xxyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_xxyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_xxyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_xxyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 239);

            auto g_x_0_xxyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_xxyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_xxyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_xxyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_xxyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_xxyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_xxyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_xxyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_xxyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_xxyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_xxyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_xxyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 251);

            #pragma omp simd aligned(cd_y, g_x_0_xxyzz_xxxxxx, g_x_0_xxyzz_xxxxxy, g_x_0_xxyzz_xxxxxz, g_x_0_xxyzz_xxxxyy, g_x_0_xxyzz_xxxxyz, g_x_0_xxyzz_xxxxzz, g_x_0_xxyzz_xxxyyy, g_x_0_xxyzz_xxxyyz, g_x_0_xxyzz_xxxyzz, g_x_0_xxyzz_xxxzzz, g_x_0_xxyzz_xxyyyy, g_x_0_xxyzz_xxyyyz, g_x_0_xxyzz_xxyyzz, g_x_0_xxyzz_xxyzzz, g_x_0_xxyzz_xxzzzz, g_x_0_xxyzz_xyyyyy, g_x_0_xxyzz_xyyyyz, g_x_0_xxyzz_xyyyzz, g_x_0_xxyzz_xyyzzz, g_x_0_xxyzz_xyzzzz, g_x_0_xxyzz_xzzzzz, g_x_0_xxyzz_yyyyyy, g_x_0_xxyzz_yyyyyz, g_x_0_xxyzz_yyyyzz, g_x_0_xxyzz_yyyzzz, g_x_0_xxyzz_yyzzzz, g_x_0_xxyzz_yzzzzz, g_x_0_xxyzz_zzzzzz, g_x_0_xxzz_xxxxxx, g_x_0_xxzz_xxxxxxy, g_x_0_xxzz_xxxxxy, g_x_0_xxzz_xxxxxyy, g_x_0_xxzz_xxxxxyz, g_x_0_xxzz_xxxxxz, g_x_0_xxzz_xxxxyy, g_x_0_xxzz_xxxxyyy, g_x_0_xxzz_xxxxyyz, g_x_0_xxzz_xxxxyz, g_x_0_xxzz_xxxxyzz, g_x_0_xxzz_xxxxzz, g_x_0_xxzz_xxxyyy, g_x_0_xxzz_xxxyyyy, g_x_0_xxzz_xxxyyyz, g_x_0_xxzz_xxxyyz, g_x_0_xxzz_xxxyyzz, g_x_0_xxzz_xxxyzz, g_x_0_xxzz_xxxyzzz, g_x_0_xxzz_xxxzzz, g_x_0_xxzz_xxyyyy, g_x_0_xxzz_xxyyyyy, g_x_0_xxzz_xxyyyyz, g_x_0_xxzz_xxyyyz, g_x_0_xxzz_xxyyyzz, g_x_0_xxzz_xxyyzz, g_x_0_xxzz_xxyyzzz, g_x_0_xxzz_xxyzzz, g_x_0_xxzz_xxyzzzz, g_x_0_xxzz_xxzzzz, g_x_0_xxzz_xyyyyy, g_x_0_xxzz_xyyyyyy, g_x_0_xxzz_xyyyyyz, g_x_0_xxzz_xyyyyz, g_x_0_xxzz_xyyyyzz, g_x_0_xxzz_xyyyzz, g_x_0_xxzz_xyyyzzz, g_x_0_xxzz_xyyzzz, g_x_0_xxzz_xyyzzzz, g_x_0_xxzz_xyzzzz, g_x_0_xxzz_xyzzzzz, g_x_0_xxzz_xzzzzz, g_x_0_xxzz_yyyyyy, g_x_0_xxzz_yyyyyyy, g_x_0_xxzz_yyyyyyz, g_x_0_xxzz_yyyyyz, g_x_0_xxzz_yyyyyzz, g_x_0_xxzz_yyyyzz, g_x_0_xxzz_yyyyzzz, g_x_0_xxzz_yyyzzz, g_x_0_xxzz_yyyzzzz, g_x_0_xxzz_yyzzzz, g_x_0_xxzz_yyzzzzz, g_x_0_xxzz_yzzzzz, g_x_0_xxzz_yzzzzzz, g_x_0_xxzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzz_xxxxxx[k] = -g_x_0_xxzz_xxxxxx[k] * cd_y[k] + g_x_0_xxzz_xxxxxxy[k];

                g_x_0_xxyzz_xxxxxy[k] = -g_x_0_xxzz_xxxxxy[k] * cd_y[k] + g_x_0_xxzz_xxxxxyy[k];

                g_x_0_xxyzz_xxxxxz[k] = -g_x_0_xxzz_xxxxxz[k] * cd_y[k] + g_x_0_xxzz_xxxxxyz[k];

                g_x_0_xxyzz_xxxxyy[k] = -g_x_0_xxzz_xxxxyy[k] * cd_y[k] + g_x_0_xxzz_xxxxyyy[k];

                g_x_0_xxyzz_xxxxyz[k] = -g_x_0_xxzz_xxxxyz[k] * cd_y[k] + g_x_0_xxzz_xxxxyyz[k];

                g_x_0_xxyzz_xxxxzz[k] = -g_x_0_xxzz_xxxxzz[k] * cd_y[k] + g_x_0_xxzz_xxxxyzz[k];

                g_x_0_xxyzz_xxxyyy[k] = -g_x_0_xxzz_xxxyyy[k] * cd_y[k] + g_x_0_xxzz_xxxyyyy[k];

                g_x_0_xxyzz_xxxyyz[k] = -g_x_0_xxzz_xxxyyz[k] * cd_y[k] + g_x_0_xxzz_xxxyyyz[k];

                g_x_0_xxyzz_xxxyzz[k] = -g_x_0_xxzz_xxxyzz[k] * cd_y[k] + g_x_0_xxzz_xxxyyzz[k];

                g_x_0_xxyzz_xxxzzz[k] = -g_x_0_xxzz_xxxzzz[k] * cd_y[k] + g_x_0_xxzz_xxxyzzz[k];

                g_x_0_xxyzz_xxyyyy[k] = -g_x_0_xxzz_xxyyyy[k] * cd_y[k] + g_x_0_xxzz_xxyyyyy[k];

                g_x_0_xxyzz_xxyyyz[k] = -g_x_0_xxzz_xxyyyz[k] * cd_y[k] + g_x_0_xxzz_xxyyyyz[k];

                g_x_0_xxyzz_xxyyzz[k] = -g_x_0_xxzz_xxyyzz[k] * cd_y[k] + g_x_0_xxzz_xxyyyzz[k];

                g_x_0_xxyzz_xxyzzz[k] = -g_x_0_xxzz_xxyzzz[k] * cd_y[k] + g_x_0_xxzz_xxyyzzz[k];

                g_x_0_xxyzz_xxzzzz[k] = -g_x_0_xxzz_xxzzzz[k] * cd_y[k] + g_x_0_xxzz_xxyzzzz[k];

                g_x_0_xxyzz_xyyyyy[k] = -g_x_0_xxzz_xyyyyy[k] * cd_y[k] + g_x_0_xxzz_xyyyyyy[k];

                g_x_0_xxyzz_xyyyyz[k] = -g_x_0_xxzz_xyyyyz[k] * cd_y[k] + g_x_0_xxzz_xyyyyyz[k];

                g_x_0_xxyzz_xyyyzz[k] = -g_x_0_xxzz_xyyyzz[k] * cd_y[k] + g_x_0_xxzz_xyyyyzz[k];

                g_x_0_xxyzz_xyyzzz[k] = -g_x_0_xxzz_xyyzzz[k] * cd_y[k] + g_x_0_xxzz_xyyyzzz[k];

                g_x_0_xxyzz_xyzzzz[k] = -g_x_0_xxzz_xyzzzz[k] * cd_y[k] + g_x_0_xxzz_xyyzzzz[k];

                g_x_0_xxyzz_xzzzzz[k] = -g_x_0_xxzz_xzzzzz[k] * cd_y[k] + g_x_0_xxzz_xyzzzzz[k];

                g_x_0_xxyzz_yyyyyy[k] = -g_x_0_xxzz_yyyyyy[k] * cd_y[k] + g_x_0_xxzz_yyyyyyy[k];

                g_x_0_xxyzz_yyyyyz[k] = -g_x_0_xxzz_yyyyyz[k] * cd_y[k] + g_x_0_xxzz_yyyyyyz[k];

                g_x_0_xxyzz_yyyyzz[k] = -g_x_0_xxzz_yyyyzz[k] * cd_y[k] + g_x_0_xxzz_yyyyyzz[k];

                g_x_0_xxyzz_yyyzzz[k] = -g_x_0_xxzz_yyyzzz[k] * cd_y[k] + g_x_0_xxzz_yyyyzzz[k];

                g_x_0_xxyzz_yyzzzz[k] = -g_x_0_xxzz_yyzzzz[k] * cd_y[k] + g_x_0_xxzz_yyyzzzz[k];

                g_x_0_xxyzz_yzzzzz[k] = -g_x_0_xxzz_yzzzzz[k] * cd_y[k] + g_x_0_xxzz_yyzzzzz[k];

                g_x_0_xxyzz_zzzzzz[k] = -g_x_0_xxzz_zzzzzz[k] * cd_y[k] + g_x_0_xxzz_yzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 252);

            auto g_x_0_xxzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_xxzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 254);

            auto g_x_0_xxzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_xxzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_xxzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 257);

            auto g_x_0_xxzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_xxzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 259);

            auto g_x_0_xxzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_xxzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 261);

            auto g_x_0_xxzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_xxzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_xxzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_xxzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_xxzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 266);

            auto g_x_0_xxzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_xxzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_xxzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 269);

            auto g_x_0_xxzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_xxzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_xxzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 272);

            auto g_x_0_xxzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 273);

            auto g_x_0_xxzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_xxzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 275);

            auto g_x_0_xxzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_xxzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_xxzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 278);

            auto g_x_0_xxzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 279);

            #pragma omp simd aligned(cd_z, g_x_0_xxzz_xxxxxx, g_x_0_xxzz_xxxxxxz, g_x_0_xxzz_xxxxxy, g_x_0_xxzz_xxxxxyz, g_x_0_xxzz_xxxxxz, g_x_0_xxzz_xxxxxzz, g_x_0_xxzz_xxxxyy, g_x_0_xxzz_xxxxyyz, g_x_0_xxzz_xxxxyz, g_x_0_xxzz_xxxxyzz, g_x_0_xxzz_xxxxzz, g_x_0_xxzz_xxxxzzz, g_x_0_xxzz_xxxyyy, g_x_0_xxzz_xxxyyyz, g_x_0_xxzz_xxxyyz, g_x_0_xxzz_xxxyyzz, g_x_0_xxzz_xxxyzz, g_x_0_xxzz_xxxyzzz, g_x_0_xxzz_xxxzzz, g_x_0_xxzz_xxxzzzz, g_x_0_xxzz_xxyyyy, g_x_0_xxzz_xxyyyyz, g_x_0_xxzz_xxyyyz, g_x_0_xxzz_xxyyyzz, g_x_0_xxzz_xxyyzz, g_x_0_xxzz_xxyyzzz, g_x_0_xxzz_xxyzzz, g_x_0_xxzz_xxyzzzz, g_x_0_xxzz_xxzzzz, g_x_0_xxzz_xxzzzzz, g_x_0_xxzz_xyyyyy, g_x_0_xxzz_xyyyyyz, g_x_0_xxzz_xyyyyz, g_x_0_xxzz_xyyyyzz, g_x_0_xxzz_xyyyzz, g_x_0_xxzz_xyyyzzz, g_x_0_xxzz_xyyzzz, g_x_0_xxzz_xyyzzzz, g_x_0_xxzz_xyzzzz, g_x_0_xxzz_xyzzzzz, g_x_0_xxzz_xzzzzz, g_x_0_xxzz_xzzzzzz, g_x_0_xxzz_yyyyyy, g_x_0_xxzz_yyyyyyz, g_x_0_xxzz_yyyyyz, g_x_0_xxzz_yyyyyzz, g_x_0_xxzz_yyyyzz, g_x_0_xxzz_yyyyzzz, g_x_0_xxzz_yyyzzz, g_x_0_xxzz_yyyzzzz, g_x_0_xxzz_yyzzzz, g_x_0_xxzz_yyzzzzz, g_x_0_xxzz_yzzzzz, g_x_0_xxzz_yzzzzzz, g_x_0_xxzz_zzzzzz, g_x_0_xxzz_zzzzzzz, g_x_0_xxzzz_xxxxxx, g_x_0_xxzzz_xxxxxy, g_x_0_xxzzz_xxxxxz, g_x_0_xxzzz_xxxxyy, g_x_0_xxzzz_xxxxyz, g_x_0_xxzzz_xxxxzz, g_x_0_xxzzz_xxxyyy, g_x_0_xxzzz_xxxyyz, g_x_0_xxzzz_xxxyzz, g_x_0_xxzzz_xxxzzz, g_x_0_xxzzz_xxyyyy, g_x_0_xxzzz_xxyyyz, g_x_0_xxzzz_xxyyzz, g_x_0_xxzzz_xxyzzz, g_x_0_xxzzz_xxzzzz, g_x_0_xxzzz_xyyyyy, g_x_0_xxzzz_xyyyyz, g_x_0_xxzzz_xyyyzz, g_x_0_xxzzz_xyyzzz, g_x_0_xxzzz_xyzzzz, g_x_0_xxzzz_xzzzzz, g_x_0_xxzzz_yyyyyy, g_x_0_xxzzz_yyyyyz, g_x_0_xxzzz_yyyyzz, g_x_0_xxzzz_yyyzzz, g_x_0_xxzzz_yyzzzz, g_x_0_xxzzz_yzzzzz, g_x_0_xxzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzz_xxxxxx[k] = -g_x_0_xxzz_xxxxxx[k] * cd_z[k] + g_x_0_xxzz_xxxxxxz[k];

                g_x_0_xxzzz_xxxxxy[k] = -g_x_0_xxzz_xxxxxy[k] * cd_z[k] + g_x_0_xxzz_xxxxxyz[k];

                g_x_0_xxzzz_xxxxxz[k] = -g_x_0_xxzz_xxxxxz[k] * cd_z[k] + g_x_0_xxzz_xxxxxzz[k];

                g_x_0_xxzzz_xxxxyy[k] = -g_x_0_xxzz_xxxxyy[k] * cd_z[k] + g_x_0_xxzz_xxxxyyz[k];

                g_x_0_xxzzz_xxxxyz[k] = -g_x_0_xxzz_xxxxyz[k] * cd_z[k] + g_x_0_xxzz_xxxxyzz[k];

                g_x_0_xxzzz_xxxxzz[k] = -g_x_0_xxzz_xxxxzz[k] * cd_z[k] + g_x_0_xxzz_xxxxzzz[k];

                g_x_0_xxzzz_xxxyyy[k] = -g_x_0_xxzz_xxxyyy[k] * cd_z[k] + g_x_0_xxzz_xxxyyyz[k];

                g_x_0_xxzzz_xxxyyz[k] = -g_x_0_xxzz_xxxyyz[k] * cd_z[k] + g_x_0_xxzz_xxxyyzz[k];

                g_x_0_xxzzz_xxxyzz[k] = -g_x_0_xxzz_xxxyzz[k] * cd_z[k] + g_x_0_xxzz_xxxyzzz[k];

                g_x_0_xxzzz_xxxzzz[k] = -g_x_0_xxzz_xxxzzz[k] * cd_z[k] + g_x_0_xxzz_xxxzzzz[k];

                g_x_0_xxzzz_xxyyyy[k] = -g_x_0_xxzz_xxyyyy[k] * cd_z[k] + g_x_0_xxzz_xxyyyyz[k];

                g_x_0_xxzzz_xxyyyz[k] = -g_x_0_xxzz_xxyyyz[k] * cd_z[k] + g_x_0_xxzz_xxyyyzz[k];

                g_x_0_xxzzz_xxyyzz[k] = -g_x_0_xxzz_xxyyzz[k] * cd_z[k] + g_x_0_xxzz_xxyyzzz[k];

                g_x_0_xxzzz_xxyzzz[k] = -g_x_0_xxzz_xxyzzz[k] * cd_z[k] + g_x_0_xxzz_xxyzzzz[k];

                g_x_0_xxzzz_xxzzzz[k] = -g_x_0_xxzz_xxzzzz[k] * cd_z[k] + g_x_0_xxzz_xxzzzzz[k];

                g_x_0_xxzzz_xyyyyy[k] = -g_x_0_xxzz_xyyyyy[k] * cd_z[k] + g_x_0_xxzz_xyyyyyz[k];

                g_x_0_xxzzz_xyyyyz[k] = -g_x_0_xxzz_xyyyyz[k] * cd_z[k] + g_x_0_xxzz_xyyyyzz[k];

                g_x_0_xxzzz_xyyyzz[k] = -g_x_0_xxzz_xyyyzz[k] * cd_z[k] + g_x_0_xxzz_xyyyzzz[k];

                g_x_0_xxzzz_xyyzzz[k] = -g_x_0_xxzz_xyyzzz[k] * cd_z[k] + g_x_0_xxzz_xyyzzzz[k];

                g_x_0_xxzzz_xyzzzz[k] = -g_x_0_xxzz_xyzzzz[k] * cd_z[k] + g_x_0_xxzz_xyzzzzz[k];

                g_x_0_xxzzz_xzzzzz[k] = -g_x_0_xxzz_xzzzzz[k] * cd_z[k] + g_x_0_xxzz_xzzzzzz[k];

                g_x_0_xxzzz_yyyyyy[k] = -g_x_0_xxzz_yyyyyy[k] * cd_z[k] + g_x_0_xxzz_yyyyyyz[k];

                g_x_0_xxzzz_yyyyyz[k] = -g_x_0_xxzz_yyyyyz[k] * cd_z[k] + g_x_0_xxzz_yyyyyzz[k];

                g_x_0_xxzzz_yyyyzz[k] = -g_x_0_xxzz_yyyyzz[k] * cd_z[k] + g_x_0_xxzz_yyyyzzz[k];

                g_x_0_xxzzz_yyyzzz[k] = -g_x_0_xxzz_yyyzzz[k] * cd_z[k] + g_x_0_xxzz_yyyzzzz[k];

                g_x_0_xxzzz_yyzzzz[k] = -g_x_0_xxzz_yyzzzz[k] * cd_z[k] + g_x_0_xxzz_yyzzzzz[k];

                g_x_0_xxzzz_yzzzzz[k] = -g_x_0_xxzz_yzzzzz[k] * cd_z[k] + g_x_0_xxzz_yzzzzzz[k];

                g_x_0_xxzzz_zzzzzz[k] = -g_x_0_xxzz_zzzzzz[k] * cd_z[k] + g_x_0_xxzz_zzzzzzz[k];
            }

            /// Set up 280-308 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 280);

            auto g_x_0_xyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 281);

            auto g_x_0_xyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 282);

            auto g_x_0_xyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 283);

            auto g_x_0_xyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 284);

            auto g_x_0_xyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 285);

            auto g_x_0_xyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 286);

            auto g_x_0_xyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 287);

            auto g_x_0_xyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 288);

            auto g_x_0_xyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 289);

            auto g_x_0_xyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 290);

            auto g_x_0_xyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 291);

            auto g_x_0_xyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 292);

            auto g_x_0_xyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 293);

            auto g_x_0_xyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 294);

            auto g_x_0_xyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 295);

            auto g_x_0_xyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 296);

            auto g_x_0_xyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 297);

            auto g_x_0_xyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 298);

            auto g_x_0_xyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 299);

            auto g_x_0_xyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 300);

            auto g_x_0_xyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 301);

            auto g_x_0_xyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 302);

            auto g_x_0_xyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 303);

            auto g_x_0_xyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 304);

            auto g_x_0_xyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 305);

            auto g_x_0_xyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 306);

            auto g_x_0_xyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 307);

            #pragma omp simd aligned(cd_y, g_x_0_xyyy_xxxxxx, g_x_0_xyyy_xxxxxxy, g_x_0_xyyy_xxxxxy, g_x_0_xyyy_xxxxxyy, g_x_0_xyyy_xxxxxyz, g_x_0_xyyy_xxxxxz, g_x_0_xyyy_xxxxyy, g_x_0_xyyy_xxxxyyy, g_x_0_xyyy_xxxxyyz, g_x_0_xyyy_xxxxyz, g_x_0_xyyy_xxxxyzz, g_x_0_xyyy_xxxxzz, g_x_0_xyyy_xxxyyy, g_x_0_xyyy_xxxyyyy, g_x_0_xyyy_xxxyyyz, g_x_0_xyyy_xxxyyz, g_x_0_xyyy_xxxyyzz, g_x_0_xyyy_xxxyzz, g_x_0_xyyy_xxxyzzz, g_x_0_xyyy_xxxzzz, g_x_0_xyyy_xxyyyy, g_x_0_xyyy_xxyyyyy, g_x_0_xyyy_xxyyyyz, g_x_0_xyyy_xxyyyz, g_x_0_xyyy_xxyyyzz, g_x_0_xyyy_xxyyzz, g_x_0_xyyy_xxyyzzz, g_x_0_xyyy_xxyzzz, g_x_0_xyyy_xxyzzzz, g_x_0_xyyy_xxzzzz, g_x_0_xyyy_xyyyyy, g_x_0_xyyy_xyyyyyy, g_x_0_xyyy_xyyyyyz, g_x_0_xyyy_xyyyyz, g_x_0_xyyy_xyyyyzz, g_x_0_xyyy_xyyyzz, g_x_0_xyyy_xyyyzzz, g_x_0_xyyy_xyyzzz, g_x_0_xyyy_xyyzzzz, g_x_0_xyyy_xyzzzz, g_x_0_xyyy_xyzzzzz, g_x_0_xyyy_xzzzzz, g_x_0_xyyy_yyyyyy, g_x_0_xyyy_yyyyyyy, g_x_0_xyyy_yyyyyyz, g_x_0_xyyy_yyyyyz, g_x_0_xyyy_yyyyyzz, g_x_0_xyyy_yyyyzz, g_x_0_xyyy_yyyyzzz, g_x_0_xyyy_yyyzzz, g_x_0_xyyy_yyyzzzz, g_x_0_xyyy_yyzzzz, g_x_0_xyyy_yyzzzzz, g_x_0_xyyy_yzzzzz, g_x_0_xyyy_yzzzzzz, g_x_0_xyyy_zzzzzz, g_x_0_xyyyy_xxxxxx, g_x_0_xyyyy_xxxxxy, g_x_0_xyyyy_xxxxxz, g_x_0_xyyyy_xxxxyy, g_x_0_xyyyy_xxxxyz, g_x_0_xyyyy_xxxxzz, g_x_0_xyyyy_xxxyyy, g_x_0_xyyyy_xxxyyz, g_x_0_xyyyy_xxxyzz, g_x_0_xyyyy_xxxzzz, g_x_0_xyyyy_xxyyyy, g_x_0_xyyyy_xxyyyz, g_x_0_xyyyy_xxyyzz, g_x_0_xyyyy_xxyzzz, g_x_0_xyyyy_xxzzzz, g_x_0_xyyyy_xyyyyy, g_x_0_xyyyy_xyyyyz, g_x_0_xyyyy_xyyyzz, g_x_0_xyyyy_xyyzzz, g_x_0_xyyyy_xyzzzz, g_x_0_xyyyy_xzzzzz, g_x_0_xyyyy_yyyyyy, g_x_0_xyyyy_yyyyyz, g_x_0_xyyyy_yyyyzz, g_x_0_xyyyy_yyyzzz, g_x_0_xyyyy_yyzzzz, g_x_0_xyyyy_yzzzzz, g_x_0_xyyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyy_xxxxxx[k] = -g_x_0_xyyy_xxxxxx[k] * cd_y[k] + g_x_0_xyyy_xxxxxxy[k];

                g_x_0_xyyyy_xxxxxy[k] = -g_x_0_xyyy_xxxxxy[k] * cd_y[k] + g_x_0_xyyy_xxxxxyy[k];

                g_x_0_xyyyy_xxxxxz[k] = -g_x_0_xyyy_xxxxxz[k] * cd_y[k] + g_x_0_xyyy_xxxxxyz[k];

                g_x_0_xyyyy_xxxxyy[k] = -g_x_0_xyyy_xxxxyy[k] * cd_y[k] + g_x_0_xyyy_xxxxyyy[k];

                g_x_0_xyyyy_xxxxyz[k] = -g_x_0_xyyy_xxxxyz[k] * cd_y[k] + g_x_0_xyyy_xxxxyyz[k];

                g_x_0_xyyyy_xxxxzz[k] = -g_x_0_xyyy_xxxxzz[k] * cd_y[k] + g_x_0_xyyy_xxxxyzz[k];

                g_x_0_xyyyy_xxxyyy[k] = -g_x_0_xyyy_xxxyyy[k] * cd_y[k] + g_x_0_xyyy_xxxyyyy[k];

                g_x_0_xyyyy_xxxyyz[k] = -g_x_0_xyyy_xxxyyz[k] * cd_y[k] + g_x_0_xyyy_xxxyyyz[k];

                g_x_0_xyyyy_xxxyzz[k] = -g_x_0_xyyy_xxxyzz[k] * cd_y[k] + g_x_0_xyyy_xxxyyzz[k];

                g_x_0_xyyyy_xxxzzz[k] = -g_x_0_xyyy_xxxzzz[k] * cd_y[k] + g_x_0_xyyy_xxxyzzz[k];

                g_x_0_xyyyy_xxyyyy[k] = -g_x_0_xyyy_xxyyyy[k] * cd_y[k] + g_x_0_xyyy_xxyyyyy[k];

                g_x_0_xyyyy_xxyyyz[k] = -g_x_0_xyyy_xxyyyz[k] * cd_y[k] + g_x_0_xyyy_xxyyyyz[k];

                g_x_0_xyyyy_xxyyzz[k] = -g_x_0_xyyy_xxyyzz[k] * cd_y[k] + g_x_0_xyyy_xxyyyzz[k];

                g_x_0_xyyyy_xxyzzz[k] = -g_x_0_xyyy_xxyzzz[k] * cd_y[k] + g_x_0_xyyy_xxyyzzz[k];

                g_x_0_xyyyy_xxzzzz[k] = -g_x_0_xyyy_xxzzzz[k] * cd_y[k] + g_x_0_xyyy_xxyzzzz[k];

                g_x_0_xyyyy_xyyyyy[k] = -g_x_0_xyyy_xyyyyy[k] * cd_y[k] + g_x_0_xyyy_xyyyyyy[k];

                g_x_0_xyyyy_xyyyyz[k] = -g_x_0_xyyy_xyyyyz[k] * cd_y[k] + g_x_0_xyyy_xyyyyyz[k];

                g_x_0_xyyyy_xyyyzz[k] = -g_x_0_xyyy_xyyyzz[k] * cd_y[k] + g_x_0_xyyy_xyyyyzz[k];

                g_x_0_xyyyy_xyyzzz[k] = -g_x_0_xyyy_xyyzzz[k] * cd_y[k] + g_x_0_xyyy_xyyyzzz[k];

                g_x_0_xyyyy_xyzzzz[k] = -g_x_0_xyyy_xyzzzz[k] * cd_y[k] + g_x_0_xyyy_xyyzzzz[k];

                g_x_0_xyyyy_xzzzzz[k] = -g_x_0_xyyy_xzzzzz[k] * cd_y[k] + g_x_0_xyyy_xyzzzzz[k];

                g_x_0_xyyyy_yyyyyy[k] = -g_x_0_xyyy_yyyyyy[k] * cd_y[k] + g_x_0_xyyy_yyyyyyy[k];

                g_x_0_xyyyy_yyyyyz[k] = -g_x_0_xyyy_yyyyyz[k] * cd_y[k] + g_x_0_xyyy_yyyyyyz[k];

                g_x_0_xyyyy_yyyyzz[k] = -g_x_0_xyyy_yyyyzz[k] * cd_y[k] + g_x_0_xyyy_yyyyyzz[k];

                g_x_0_xyyyy_yyyzzz[k] = -g_x_0_xyyy_yyyzzz[k] * cd_y[k] + g_x_0_xyyy_yyyyzzz[k];

                g_x_0_xyyyy_yyzzzz[k] = -g_x_0_xyyy_yyzzzz[k] * cd_y[k] + g_x_0_xyyy_yyyzzzz[k];

                g_x_0_xyyyy_yzzzzz[k] = -g_x_0_xyyy_yzzzzz[k] * cd_y[k] + g_x_0_xyyy_yyzzzzz[k];

                g_x_0_xyyyy_zzzzzz[k] = -g_x_0_xyyy_zzzzzz[k] * cd_y[k] + g_x_0_xyyy_yzzzzzz[k];
            }

            /// Set up 308-336 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 308);

            auto g_x_0_xyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 309);

            auto g_x_0_xyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 310);

            auto g_x_0_xyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 311);

            auto g_x_0_xyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 312);

            auto g_x_0_xyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 313);

            auto g_x_0_xyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 314);

            auto g_x_0_xyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 315);

            auto g_x_0_xyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 316);

            auto g_x_0_xyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 317);

            auto g_x_0_xyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 318);

            auto g_x_0_xyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 319);

            auto g_x_0_xyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 320);

            auto g_x_0_xyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 321);

            auto g_x_0_xyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 322);

            auto g_x_0_xyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 323);

            auto g_x_0_xyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 324);

            auto g_x_0_xyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 325);

            auto g_x_0_xyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 326);

            auto g_x_0_xyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 327);

            auto g_x_0_xyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 328);

            auto g_x_0_xyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 329);

            auto g_x_0_xyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 330);

            auto g_x_0_xyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 331);

            auto g_x_0_xyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 332);

            auto g_x_0_xyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 333);

            auto g_x_0_xyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 334);

            auto g_x_0_xyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 335);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyz_xxxxxx, g_x_0_xyyyz_xxxxxy, g_x_0_xyyyz_xxxxxz, g_x_0_xyyyz_xxxxyy, g_x_0_xyyyz_xxxxyz, g_x_0_xyyyz_xxxxzz, g_x_0_xyyyz_xxxyyy, g_x_0_xyyyz_xxxyyz, g_x_0_xyyyz_xxxyzz, g_x_0_xyyyz_xxxzzz, g_x_0_xyyyz_xxyyyy, g_x_0_xyyyz_xxyyyz, g_x_0_xyyyz_xxyyzz, g_x_0_xyyyz_xxyzzz, g_x_0_xyyyz_xxzzzz, g_x_0_xyyyz_xyyyyy, g_x_0_xyyyz_xyyyyz, g_x_0_xyyyz_xyyyzz, g_x_0_xyyyz_xyyzzz, g_x_0_xyyyz_xyzzzz, g_x_0_xyyyz_xzzzzz, g_x_0_xyyyz_yyyyyy, g_x_0_xyyyz_yyyyyz, g_x_0_xyyyz_yyyyzz, g_x_0_xyyyz_yyyzzz, g_x_0_xyyyz_yyzzzz, g_x_0_xyyyz_yzzzzz, g_x_0_xyyyz_zzzzzz, g_x_0_xyyz_xxxxxx, g_x_0_xyyz_xxxxxxy, g_x_0_xyyz_xxxxxy, g_x_0_xyyz_xxxxxyy, g_x_0_xyyz_xxxxxyz, g_x_0_xyyz_xxxxxz, g_x_0_xyyz_xxxxyy, g_x_0_xyyz_xxxxyyy, g_x_0_xyyz_xxxxyyz, g_x_0_xyyz_xxxxyz, g_x_0_xyyz_xxxxyzz, g_x_0_xyyz_xxxxzz, g_x_0_xyyz_xxxyyy, g_x_0_xyyz_xxxyyyy, g_x_0_xyyz_xxxyyyz, g_x_0_xyyz_xxxyyz, g_x_0_xyyz_xxxyyzz, g_x_0_xyyz_xxxyzz, g_x_0_xyyz_xxxyzzz, g_x_0_xyyz_xxxzzz, g_x_0_xyyz_xxyyyy, g_x_0_xyyz_xxyyyyy, g_x_0_xyyz_xxyyyyz, g_x_0_xyyz_xxyyyz, g_x_0_xyyz_xxyyyzz, g_x_0_xyyz_xxyyzz, g_x_0_xyyz_xxyyzzz, g_x_0_xyyz_xxyzzz, g_x_0_xyyz_xxyzzzz, g_x_0_xyyz_xxzzzz, g_x_0_xyyz_xyyyyy, g_x_0_xyyz_xyyyyyy, g_x_0_xyyz_xyyyyyz, g_x_0_xyyz_xyyyyz, g_x_0_xyyz_xyyyyzz, g_x_0_xyyz_xyyyzz, g_x_0_xyyz_xyyyzzz, g_x_0_xyyz_xyyzzz, g_x_0_xyyz_xyyzzzz, g_x_0_xyyz_xyzzzz, g_x_0_xyyz_xyzzzzz, g_x_0_xyyz_xzzzzz, g_x_0_xyyz_yyyyyy, g_x_0_xyyz_yyyyyyy, g_x_0_xyyz_yyyyyyz, g_x_0_xyyz_yyyyyz, g_x_0_xyyz_yyyyyzz, g_x_0_xyyz_yyyyzz, g_x_0_xyyz_yyyyzzz, g_x_0_xyyz_yyyzzz, g_x_0_xyyz_yyyzzzz, g_x_0_xyyz_yyzzzz, g_x_0_xyyz_yyzzzzz, g_x_0_xyyz_yzzzzz, g_x_0_xyyz_yzzzzzz, g_x_0_xyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyz_xxxxxx[k] = -g_x_0_xyyz_xxxxxx[k] * cd_y[k] + g_x_0_xyyz_xxxxxxy[k];

                g_x_0_xyyyz_xxxxxy[k] = -g_x_0_xyyz_xxxxxy[k] * cd_y[k] + g_x_0_xyyz_xxxxxyy[k];

                g_x_0_xyyyz_xxxxxz[k] = -g_x_0_xyyz_xxxxxz[k] * cd_y[k] + g_x_0_xyyz_xxxxxyz[k];

                g_x_0_xyyyz_xxxxyy[k] = -g_x_0_xyyz_xxxxyy[k] * cd_y[k] + g_x_0_xyyz_xxxxyyy[k];

                g_x_0_xyyyz_xxxxyz[k] = -g_x_0_xyyz_xxxxyz[k] * cd_y[k] + g_x_0_xyyz_xxxxyyz[k];

                g_x_0_xyyyz_xxxxzz[k] = -g_x_0_xyyz_xxxxzz[k] * cd_y[k] + g_x_0_xyyz_xxxxyzz[k];

                g_x_0_xyyyz_xxxyyy[k] = -g_x_0_xyyz_xxxyyy[k] * cd_y[k] + g_x_0_xyyz_xxxyyyy[k];

                g_x_0_xyyyz_xxxyyz[k] = -g_x_0_xyyz_xxxyyz[k] * cd_y[k] + g_x_0_xyyz_xxxyyyz[k];

                g_x_0_xyyyz_xxxyzz[k] = -g_x_0_xyyz_xxxyzz[k] * cd_y[k] + g_x_0_xyyz_xxxyyzz[k];

                g_x_0_xyyyz_xxxzzz[k] = -g_x_0_xyyz_xxxzzz[k] * cd_y[k] + g_x_0_xyyz_xxxyzzz[k];

                g_x_0_xyyyz_xxyyyy[k] = -g_x_0_xyyz_xxyyyy[k] * cd_y[k] + g_x_0_xyyz_xxyyyyy[k];

                g_x_0_xyyyz_xxyyyz[k] = -g_x_0_xyyz_xxyyyz[k] * cd_y[k] + g_x_0_xyyz_xxyyyyz[k];

                g_x_0_xyyyz_xxyyzz[k] = -g_x_0_xyyz_xxyyzz[k] * cd_y[k] + g_x_0_xyyz_xxyyyzz[k];

                g_x_0_xyyyz_xxyzzz[k] = -g_x_0_xyyz_xxyzzz[k] * cd_y[k] + g_x_0_xyyz_xxyyzzz[k];

                g_x_0_xyyyz_xxzzzz[k] = -g_x_0_xyyz_xxzzzz[k] * cd_y[k] + g_x_0_xyyz_xxyzzzz[k];

                g_x_0_xyyyz_xyyyyy[k] = -g_x_0_xyyz_xyyyyy[k] * cd_y[k] + g_x_0_xyyz_xyyyyyy[k];

                g_x_0_xyyyz_xyyyyz[k] = -g_x_0_xyyz_xyyyyz[k] * cd_y[k] + g_x_0_xyyz_xyyyyyz[k];

                g_x_0_xyyyz_xyyyzz[k] = -g_x_0_xyyz_xyyyzz[k] * cd_y[k] + g_x_0_xyyz_xyyyyzz[k];

                g_x_0_xyyyz_xyyzzz[k] = -g_x_0_xyyz_xyyzzz[k] * cd_y[k] + g_x_0_xyyz_xyyyzzz[k];

                g_x_0_xyyyz_xyzzzz[k] = -g_x_0_xyyz_xyzzzz[k] * cd_y[k] + g_x_0_xyyz_xyyzzzz[k];

                g_x_0_xyyyz_xzzzzz[k] = -g_x_0_xyyz_xzzzzz[k] * cd_y[k] + g_x_0_xyyz_xyzzzzz[k];

                g_x_0_xyyyz_yyyyyy[k] = -g_x_0_xyyz_yyyyyy[k] * cd_y[k] + g_x_0_xyyz_yyyyyyy[k];

                g_x_0_xyyyz_yyyyyz[k] = -g_x_0_xyyz_yyyyyz[k] * cd_y[k] + g_x_0_xyyz_yyyyyyz[k];

                g_x_0_xyyyz_yyyyzz[k] = -g_x_0_xyyz_yyyyzz[k] * cd_y[k] + g_x_0_xyyz_yyyyyzz[k];

                g_x_0_xyyyz_yyyzzz[k] = -g_x_0_xyyz_yyyzzz[k] * cd_y[k] + g_x_0_xyyz_yyyyzzz[k];

                g_x_0_xyyyz_yyzzzz[k] = -g_x_0_xyyz_yyzzzz[k] * cd_y[k] + g_x_0_xyyz_yyyzzzz[k];

                g_x_0_xyyyz_yzzzzz[k] = -g_x_0_xyyz_yzzzzz[k] * cd_y[k] + g_x_0_xyyz_yyzzzzz[k];

                g_x_0_xyyyz_zzzzzz[k] = -g_x_0_xyyz_zzzzzz[k] * cd_y[k] + g_x_0_xyyz_yzzzzzz[k];
            }

            /// Set up 336-364 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 336);

            auto g_x_0_xyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 337);

            auto g_x_0_xyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 338);

            auto g_x_0_xyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 339);

            auto g_x_0_xyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 340);

            auto g_x_0_xyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 341);

            auto g_x_0_xyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 342);

            auto g_x_0_xyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 343);

            auto g_x_0_xyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 344);

            auto g_x_0_xyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 345);

            auto g_x_0_xyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 346);

            auto g_x_0_xyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 347);

            auto g_x_0_xyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 348);

            auto g_x_0_xyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 349);

            auto g_x_0_xyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 350);

            auto g_x_0_xyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 351);

            auto g_x_0_xyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 352);

            auto g_x_0_xyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 353);

            auto g_x_0_xyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 354);

            auto g_x_0_xyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 355);

            auto g_x_0_xyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 356);

            auto g_x_0_xyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 357);

            auto g_x_0_xyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 358);

            auto g_x_0_xyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 359);

            auto g_x_0_xyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 360);

            auto g_x_0_xyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 361);

            auto g_x_0_xyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 362);

            auto g_x_0_xyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 363);

            #pragma omp simd aligned(cd_y, g_x_0_xyyzz_xxxxxx, g_x_0_xyyzz_xxxxxy, g_x_0_xyyzz_xxxxxz, g_x_0_xyyzz_xxxxyy, g_x_0_xyyzz_xxxxyz, g_x_0_xyyzz_xxxxzz, g_x_0_xyyzz_xxxyyy, g_x_0_xyyzz_xxxyyz, g_x_0_xyyzz_xxxyzz, g_x_0_xyyzz_xxxzzz, g_x_0_xyyzz_xxyyyy, g_x_0_xyyzz_xxyyyz, g_x_0_xyyzz_xxyyzz, g_x_0_xyyzz_xxyzzz, g_x_0_xyyzz_xxzzzz, g_x_0_xyyzz_xyyyyy, g_x_0_xyyzz_xyyyyz, g_x_0_xyyzz_xyyyzz, g_x_0_xyyzz_xyyzzz, g_x_0_xyyzz_xyzzzz, g_x_0_xyyzz_xzzzzz, g_x_0_xyyzz_yyyyyy, g_x_0_xyyzz_yyyyyz, g_x_0_xyyzz_yyyyzz, g_x_0_xyyzz_yyyzzz, g_x_0_xyyzz_yyzzzz, g_x_0_xyyzz_yzzzzz, g_x_0_xyyzz_zzzzzz, g_x_0_xyzz_xxxxxx, g_x_0_xyzz_xxxxxxy, g_x_0_xyzz_xxxxxy, g_x_0_xyzz_xxxxxyy, g_x_0_xyzz_xxxxxyz, g_x_0_xyzz_xxxxxz, g_x_0_xyzz_xxxxyy, g_x_0_xyzz_xxxxyyy, g_x_0_xyzz_xxxxyyz, g_x_0_xyzz_xxxxyz, g_x_0_xyzz_xxxxyzz, g_x_0_xyzz_xxxxzz, g_x_0_xyzz_xxxyyy, g_x_0_xyzz_xxxyyyy, g_x_0_xyzz_xxxyyyz, g_x_0_xyzz_xxxyyz, g_x_0_xyzz_xxxyyzz, g_x_0_xyzz_xxxyzz, g_x_0_xyzz_xxxyzzz, g_x_0_xyzz_xxxzzz, g_x_0_xyzz_xxyyyy, g_x_0_xyzz_xxyyyyy, g_x_0_xyzz_xxyyyyz, g_x_0_xyzz_xxyyyz, g_x_0_xyzz_xxyyyzz, g_x_0_xyzz_xxyyzz, g_x_0_xyzz_xxyyzzz, g_x_0_xyzz_xxyzzz, g_x_0_xyzz_xxyzzzz, g_x_0_xyzz_xxzzzz, g_x_0_xyzz_xyyyyy, g_x_0_xyzz_xyyyyyy, g_x_0_xyzz_xyyyyyz, g_x_0_xyzz_xyyyyz, g_x_0_xyzz_xyyyyzz, g_x_0_xyzz_xyyyzz, g_x_0_xyzz_xyyyzzz, g_x_0_xyzz_xyyzzz, g_x_0_xyzz_xyyzzzz, g_x_0_xyzz_xyzzzz, g_x_0_xyzz_xyzzzzz, g_x_0_xyzz_xzzzzz, g_x_0_xyzz_yyyyyy, g_x_0_xyzz_yyyyyyy, g_x_0_xyzz_yyyyyyz, g_x_0_xyzz_yyyyyz, g_x_0_xyzz_yyyyyzz, g_x_0_xyzz_yyyyzz, g_x_0_xyzz_yyyyzzz, g_x_0_xyzz_yyyzzz, g_x_0_xyzz_yyyzzzz, g_x_0_xyzz_yyzzzz, g_x_0_xyzz_yyzzzzz, g_x_0_xyzz_yzzzzz, g_x_0_xyzz_yzzzzzz, g_x_0_xyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzz_xxxxxx[k] = -g_x_0_xyzz_xxxxxx[k] * cd_y[k] + g_x_0_xyzz_xxxxxxy[k];

                g_x_0_xyyzz_xxxxxy[k] = -g_x_0_xyzz_xxxxxy[k] * cd_y[k] + g_x_0_xyzz_xxxxxyy[k];

                g_x_0_xyyzz_xxxxxz[k] = -g_x_0_xyzz_xxxxxz[k] * cd_y[k] + g_x_0_xyzz_xxxxxyz[k];

                g_x_0_xyyzz_xxxxyy[k] = -g_x_0_xyzz_xxxxyy[k] * cd_y[k] + g_x_0_xyzz_xxxxyyy[k];

                g_x_0_xyyzz_xxxxyz[k] = -g_x_0_xyzz_xxxxyz[k] * cd_y[k] + g_x_0_xyzz_xxxxyyz[k];

                g_x_0_xyyzz_xxxxzz[k] = -g_x_0_xyzz_xxxxzz[k] * cd_y[k] + g_x_0_xyzz_xxxxyzz[k];

                g_x_0_xyyzz_xxxyyy[k] = -g_x_0_xyzz_xxxyyy[k] * cd_y[k] + g_x_0_xyzz_xxxyyyy[k];

                g_x_0_xyyzz_xxxyyz[k] = -g_x_0_xyzz_xxxyyz[k] * cd_y[k] + g_x_0_xyzz_xxxyyyz[k];

                g_x_0_xyyzz_xxxyzz[k] = -g_x_0_xyzz_xxxyzz[k] * cd_y[k] + g_x_0_xyzz_xxxyyzz[k];

                g_x_0_xyyzz_xxxzzz[k] = -g_x_0_xyzz_xxxzzz[k] * cd_y[k] + g_x_0_xyzz_xxxyzzz[k];

                g_x_0_xyyzz_xxyyyy[k] = -g_x_0_xyzz_xxyyyy[k] * cd_y[k] + g_x_0_xyzz_xxyyyyy[k];

                g_x_0_xyyzz_xxyyyz[k] = -g_x_0_xyzz_xxyyyz[k] * cd_y[k] + g_x_0_xyzz_xxyyyyz[k];

                g_x_0_xyyzz_xxyyzz[k] = -g_x_0_xyzz_xxyyzz[k] * cd_y[k] + g_x_0_xyzz_xxyyyzz[k];

                g_x_0_xyyzz_xxyzzz[k] = -g_x_0_xyzz_xxyzzz[k] * cd_y[k] + g_x_0_xyzz_xxyyzzz[k];

                g_x_0_xyyzz_xxzzzz[k] = -g_x_0_xyzz_xxzzzz[k] * cd_y[k] + g_x_0_xyzz_xxyzzzz[k];

                g_x_0_xyyzz_xyyyyy[k] = -g_x_0_xyzz_xyyyyy[k] * cd_y[k] + g_x_0_xyzz_xyyyyyy[k];

                g_x_0_xyyzz_xyyyyz[k] = -g_x_0_xyzz_xyyyyz[k] * cd_y[k] + g_x_0_xyzz_xyyyyyz[k];

                g_x_0_xyyzz_xyyyzz[k] = -g_x_0_xyzz_xyyyzz[k] * cd_y[k] + g_x_0_xyzz_xyyyyzz[k];

                g_x_0_xyyzz_xyyzzz[k] = -g_x_0_xyzz_xyyzzz[k] * cd_y[k] + g_x_0_xyzz_xyyyzzz[k];

                g_x_0_xyyzz_xyzzzz[k] = -g_x_0_xyzz_xyzzzz[k] * cd_y[k] + g_x_0_xyzz_xyyzzzz[k];

                g_x_0_xyyzz_xzzzzz[k] = -g_x_0_xyzz_xzzzzz[k] * cd_y[k] + g_x_0_xyzz_xyzzzzz[k];

                g_x_0_xyyzz_yyyyyy[k] = -g_x_0_xyzz_yyyyyy[k] * cd_y[k] + g_x_0_xyzz_yyyyyyy[k];

                g_x_0_xyyzz_yyyyyz[k] = -g_x_0_xyzz_yyyyyz[k] * cd_y[k] + g_x_0_xyzz_yyyyyyz[k];

                g_x_0_xyyzz_yyyyzz[k] = -g_x_0_xyzz_yyyyzz[k] * cd_y[k] + g_x_0_xyzz_yyyyyzz[k];

                g_x_0_xyyzz_yyyzzz[k] = -g_x_0_xyzz_yyyzzz[k] * cd_y[k] + g_x_0_xyzz_yyyyzzz[k];

                g_x_0_xyyzz_yyzzzz[k] = -g_x_0_xyzz_yyzzzz[k] * cd_y[k] + g_x_0_xyzz_yyyzzzz[k];

                g_x_0_xyyzz_yzzzzz[k] = -g_x_0_xyzz_yzzzzz[k] * cd_y[k] + g_x_0_xyzz_yyzzzzz[k];

                g_x_0_xyyzz_zzzzzz[k] = -g_x_0_xyzz_zzzzzz[k] * cd_y[k] + g_x_0_xyzz_yzzzzzz[k];
            }

            /// Set up 364-392 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 364);

            auto g_x_0_xyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 365);

            auto g_x_0_xyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 366);

            auto g_x_0_xyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 367);

            auto g_x_0_xyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 368);

            auto g_x_0_xyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 369);

            auto g_x_0_xyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 370);

            auto g_x_0_xyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 371);

            auto g_x_0_xyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 372);

            auto g_x_0_xyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 373);

            auto g_x_0_xyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 374);

            auto g_x_0_xyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 375);

            auto g_x_0_xyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 376);

            auto g_x_0_xyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 377);

            auto g_x_0_xyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 378);

            auto g_x_0_xyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 379);

            auto g_x_0_xyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 380);

            auto g_x_0_xyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 381);

            auto g_x_0_xyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 382);

            auto g_x_0_xyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 383);

            auto g_x_0_xyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 384);

            auto g_x_0_xyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 385);

            auto g_x_0_xyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 386);

            auto g_x_0_xyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 387);

            auto g_x_0_xyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 388);

            auto g_x_0_xyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 389);

            auto g_x_0_xyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 390);

            auto g_x_0_xyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 391);

            #pragma omp simd aligned(cd_y, g_x_0_xyzzz_xxxxxx, g_x_0_xyzzz_xxxxxy, g_x_0_xyzzz_xxxxxz, g_x_0_xyzzz_xxxxyy, g_x_0_xyzzz_xxxxyz, g_x_0_xyzzz_xxxxzz, g_x_0_xyzzz_xxxyyy, g_x_0_xyzzz_xxxyyz, g_x_0_xyzzz_xxxyzz, g_x_0_xyzzz_xxxzzz, g_x_0_xyzzz_xxyyyy, g_x_0_xyzzz_xxyyyz, g_x_0_xyzzz_xxyyzz, g_x_0_xyzzz_xxyzzz, g_x_0_xyzzz_xxzzzz, g_x_0_xyzzz_xyyyyy, g_x_0_xyzzz_xyyyyz, g_x_0_xyzzz_xyyyzz, g_x_0_xyzzz_xyyzzz, g_x_0_xyzzz_xyzzzz, g_x_0_xyzzz_xzzzzz, g_x_0_xyzzz_yyyyyy, g_x_0_xyzzz_yyyyyz, g_x_0_xyzzz_yyyyzz, g_x_0_xyzzz_yyyzzz, g_x_0_xyzzz_yyzzzz, g_x_0_xyzzz_yzzzzz, g_x_0_xyzzz_zzzzzz, g_x_0_xzzz_xxxxxx, g_x_0_xzzz_xxxxxxy, g_x_0_xzzz_xxxxxy, g_x_0_xzzz_xxxxxyy, g_x_0_xzzz_xxxxxyz, g_x_0_xzzz_xxxxxz, g_x_0_xzzz_xxxxyy, g_x_0_xzzz_xxxxyyy, g_x_0_xzzz_xxxxyyz, g_x_0_xzzz_xxxxyz, g_x_0_xzzz_xxxxyzz, g_x_0_xzzz_xxxxzz, g_x_0_xzzz_xxxyyy, g_x_0_xzzz_xxxyyyy, g_x_0_xzzz_xxxyyyz, g_x_0_xzzz_xxxyyz, g_x_0_xzzz_xxxyyzz, g_x_0_xzzz_xxxyzz, g_x_0_xzzz_xxxyzzz, g_x_0_xzzz_xxxzzz, g_x_0_xzzz_xxyyyy, g_x_0_xzzz_xxyyyyy, g_x_0_xzzz_xxyyyyz, g_x_0_xzzz_xxyyyz, g_x_0_xzzz_xxyyyzz, g_x_0_xzzz_xxyyzz, g_x_0_xzzz_xxyyzzz, g_x_0_xzzz_xxyzzz, g_x_0_xzzz_xxyzzzz, g_x_0_xzzz_xxzzzz, g_x_0_xzzz_xyyyyy, g_x_0_xzzz_xyyyyyy, g_x_0_xzzz_xyyyyyz, g_x_0_xzzz_xyyyyz, g_x_0_xzzz_xyyyyzz, g_x_0_xzzz_xyyyzz, g_x_0_xzzz_xyyyzzz, g_x_0_xzzz_xyyzzz, g_x_0_xzzz_xyyzzzz, g_x_0_xzzz_xyzzzz, g_x_0_xzzz_xyzzzzz, g_x_0_xzzz_xzzzzz, g_x_0_xzzz_yyyyyy, g_x_0_xzzz_yyyyyyy, g_x_0_xzzz_yyyyyyz, g_x_0_xzzz_yyyyyz, g_x_0_xzzz_yyyyyzz, g_x_0_xzzz_yyyyzz, g_x_0_xzzz_yyyyzzz, g_x_0_xzzz_yyyzzz, g_x_0_xzzz_yyyzzzz, g_x_0_xzzz_yyzzzz, g_x_0_xzzz_yyzzzzz, g_x_0_xzzz_yzzzzz, g_x_0_xzzz_yzzzzzz, g_x_0_xzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzz_xxxxxx[k] = -g_x_0_xzzz_xxxxxx[k] * cd_y[k] + g_x_0_xzzz_xxxxxxy[k];

                g_x_0_xyzzz_xxxxxy[k] = -g_x_0_xzzz_xxxxxy[k] * cd_y[k] + g_x_0_xzzz_xxxxxyy[k];

                g_x_0_xyzzz_xxxxxz[k] = -g_x_0_xzzz_xxxxxz[k] * cd_y[k] + g_x_0_xzzz_xxxxxyz[k];

                g_x_0_xyzzz_xxxxyy[k] = -g_x_0_xzzz_xxxxyy[k] * cd_y[k] + g_x_0_xzzz_xxxxyyy[k];

                g_x_0_xyzzz_xxxxyz[k] = -g_x_0_xzzz_xxxxyz[k] * cd_y[k] + g_x_0_xzzz_xxxxyyz[k];

                g_x_0_xyzzz_xxxxzz[k] = -g_x_0_xzzz_xxxxzz[k] * cd_y[k] + g_x_0_xzzz_xxxxyzz[k];

                g_x_0_xyzzz_xxxyyy[k] = -g_x_0_xzzz_xxxyyy[k] * cd_y[k] + g_x_0_xzzz_xxxyyyy[k];

                g_x_0_xyzzz_xxxyyz[k] = -g_x_0_xzzz_xxxyyz[k] * cd_y[k] + g_x_0_xzzz_xxxyyyz[k];

                g_x_0_xyzzz_xxxyzz[k] = -g_x_0_xzzz_xxxyzz[k] * cd_y[k] + g_x_0_xzzz_xxxyyzz[k];

                g_x_0_xyzzz_xxxzzz[k] = -g_x_0_xzzz_xxxzzz[k] * cd_y[k] + g_x_0_xzzz_xxxyzzz[k];

                g_x_0_xyzzz_xxyyyy[k] = -g_x_0_xzzz_xxyyyy[k] * cd_y[k] + g_x_0_xzzz_xxyyyyy[k];

                g_x_0_xyzzz_xxyyyz[k] = -g_x_0_xzzz_xxyyyz[k] * cd_y[k] + g_x_0_xzzz_xxyyyyz[k];

                g_x_0_xyzzz_xxyyzz[k] = -g_x_0_xzzz_xxyyzz[k] * cd_y[k] + g_x_0_xzzz_xxyyyzz[k];

                g_x_0_xyzzz_xxyzzz[k] = -g_x_0_xzzz_xxyzzz[k] * cd_y[k] + g_x_0_xzzz_xxyyzzz[k];

                g_x_0_xyzzz_xxzzzz[k] = -g_x_0_xzzz_xxzzzz[k] * cd_y[k] + g_x_0_xzzz_xxyzzzz[k];

                g_x_0_xyzzz_xyyyyy[k] = -g_x_0_xzzz_xyyyyy[k] * cd_y[k] + g_x_0_xzzz_xyyyyyy[k];

                g_x_0_xyzzz_xyyyyz[k] = -g_x_0_xzzz_xyyyyz[k] * cd_y[k] + g_x_0_xzzz_xyyyyyz[k];

                g_x_0_xyzzz_xyyyzz[k] = -g_x_0_xzzz_xyyyzz[k] * cd_y[k] + g_x_0_xzzz_xyyyyzz[k];

                g_x_0_xyzzz_xyyzzz[k] = -g_x_0_xzzz_xyyzzz[k] * cd_y[k] + g_x_0_xzzz_xyyyzzz[k];

                g_x_0_xyzzz_xyzzzz[k] = -g_x_0_xzzz_xyzzzz[k] * cd_y[k] + g_x_0_xzzz_xyyzzzz[k];

                g_x_0_xyzzz_xzzzzz[k] = -g_x_0_xzzz_xzzzzz[k] * cd_y[k] + g_x_0_xzzz_xyzzzzz[k];

                g_x_0_xyzzz_yyyyyy[k] = -g_x_0_xzzz_yyyyyy[k] * cd_y[k] + g_x_0_xzzz_yyyyyyy[k];

                g_x_0_xyzzz_yyyyyz[k] = -g_x_0_xzzz_yyyyyz[k] * cd_y[k] + g_x_0_xzzz_yyyyyyz[k];

                g_x_0_xyzzz_yyyyzz[k] = -g_x_0_xzzz_yyyyzz[k] * cd_y[k] + g_x_0_xzzz_yyyyyzz[k];

                g_x_0_xyzzz_yyyzzz[k] = -g_x_0_xzzz_yyyzzz[k] * cd_y[k] + g_x_0_xzzz_yyyyzzz[k];

                g_x_0_xyzzz_yyzzzz[k] = -g_x_0_xzzz_yyzzzz[k] * cd_y[k] + g_x_0_xzzz_yyyzzzz[k];

                g_x_0_xyzzz_yzzzzz[k] = -g_x_0_xzzz_yzzzzz[k] * cd_y[k] + g_x_0_xzzz_yyzzzzz[k];

                g_x_0_xyzzz_zzzzzz[k] = -g_x_0_xzzz_zzzzzz[k] * cd_y[k] + g_x_0_xzzz_yzzzzzz[k];
            }

            /// Set up 392-420 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 392);

            auto g_x_0_xzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 393);

            auto g_x_0_xzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 394);

            auto g_x_0_xzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 395);

            auto g_x_0_xzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 396);

            auto g_x_0_xzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 397);

            auto g_x_0_xzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 398);

            auto g_x_0_xzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 399);

            auto g_x_0_xzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 400);

            auto g_x_0_xzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 401);

            auto g_x_0_xzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 402);

            auto g_x_0_xzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 403);

            auto g_x_0_xzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 404);

            auto g_x_0_xzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 405);

            auto g_x_0_xzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 406);

            auto g_x_0_xzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 407);

            auto g_x_0_xzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 408);

            auto g_x_0_xzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 409);

            auto g_x_0_xzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 410);

            auto g_x_0_xzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 411);

            auto g_x_0_xzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 412);

            auto g_x_0_xzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 413);

            auto g_x_0_xzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 414);

            auto g_x_0_xzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 415);

            auto g_x_0_xzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 416);

            auto g_x_0_xzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 417);

            auto g_x_0_xzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 418);

            auto g_x_0_xzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 419);

            #pragma omp simd aligned(cd_z, g_x_0_xzzz_xxxxxx, g_x_0_xzzz_xxxxxxz, g_x_0_xzzz_xxxxxy, g_x_0_xzzz_xxxxxyz, g_x_0_xzzz_xxxxxz, g_x_0_xzzz_xxxxxzz, g_x_0_xzzz_xxxxyy, g_x_0_xzzz_xxxxyyz, g_x_0_xzzz_xxxxyz, g_x_0_xzzz_xxxxyzz, g_x_0_xzzz_xxxxzz, g_x_0_xzzz_xxxxzzz, g_x_0_xzzz_xxxyyy, g_x_0_xzzz_xxxyyyz, g_x_0_xzzz_xxxyyz, g_x_0_xzzz_xxxyyzz, g_x_0_xzzz_xxxyzz, g_x_0_xzzz_xxxyzzz, g_x_0_xzzz_xxxzzz, g_x_0_xzzz_xxxzzzz, g_x_0_xzzz_xxyyyy, g_x_0_xzzz_xxyyyyz, g_x_0_xzzz_xxyyyz, g_x_0_xzzz_xxyyyzz, g_x_0_xzzz_xxyyzz, g_x_0_xzzz_xxyyzzz, g_x_0_xzzz_xxyzzz, g_x_0_xzzz_xxyzzzz, g_x_0_xzzz_xxzzzz, g_x_0_xzzz_xxzzzzz, g_x_0_xzzz_xyyyyy, g_x_0_xzzz_xyyyyyz, g_x_0_xzzz_xyyyyz, g_x_0_xzzz_xyyyyzz, g_x_0_xzzz_xyyyzz, g_x_0_xzzz_xyyyzzz, g_x_0_xzzz_xyyzzz, g_x_0_xzzz_xyyzzzz, g_x_0_xzzz_xyzzzz, g_x_0_xzzz_xyzzzzz, g_x_0_xzzz_xzzzzz, g_x_0_xzzz_xzzzzzz, g_x_0_xzzz_yyyyyy, g_x_0_xzzz_yyyyyyz, g_x_0_xzzz_yyyyyz, g_x_0_xzzz_yyyyyzz, g_x_0_xzzz_yyyyzz, g_x_0_xzzz_yyyyzzz, g_x_0_xzzz_yyyzzz, g_x_0_xzzz_yyyzzzz, g_x_0_xzzz_yyzzzz, g_x_0_xzzz_yyzzzzz, g_x_0_xzzz_yzzzzz, g_x_0_xzzz_yzzzzzz, g_x_0_xzzz_zzzzzz, g_x_0_xzzz_zzzzzzz, g_x_0_xzzzz_xxxxxx, g_x_0_xzzzz_xxxxxy, g_x_0_xzzzz_xxxxxz, g_x_0_xzzzz_xxxxyy, g_x_0_xzzzz_xxxxyz, g_x_0_xzzzz_xxxxzz, g_x_0_xzzzz_xxxyyy, g_x_0_xzzzz_xxxyyz, g_x_0_xzzzz_xxxyzz, g_x_0_xzzzz_xxxzzz, g_x_0_xzzzz_xxyyyy, g_x_0_xzzzz_xxyyyz, g_x_0_xzzzz_xxyyzz, g_x_0_xzzzz_xxyzzz, g_x_0_xzzzz_xxzzzz, g_x_0_xzzzz_xyyyyy, g_x_0_xzzzz_xyyyyz, g_x_0_xzzzz_xyyyzz, g_x_0_xzzzz_xyyzzz, g_x_0_xzzzz_xyzzzz, g_x_0_xzzzz_xzzzzz, g_x_0_xzzzz_yyyyyy, g_x_0_xzzzz_yyyyyz, g_x_0_xzzzz_yyyyzz, g_x_0_xzzzz_yyyzzz, g_x_0_xzzzz_yyzzzz, g_x_0_xzzzz_yzzzzz, g_x_0_xzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzz_xxxxxx[k] = -g_x_0_xzzz_xxxxxx[k] * cd_z[k] + g_x_0_xzzz_xxxxxxz[k];

                g_x_0_xzzzz_xxxxxy[k] = -g_x_0_xzzz_xxxxxy[k] * cd_z[k] + g_x_0_xzzz_xxxxxyz[k];

                g_x_0_xzzzz_xxxxxz[k] = -g_x_0_xzzz_xxxxxz[k] * cd_z[k] + g_x_0_xzzz_xxxxxzz[k];

                g_x_0_xzzzz_xxxxyy[k] = -g_x_0_xzzz_xxxxyy[k] * cd_z[k] + g_x_0_xzzz_xxxxyyz[k];

                g_x_0_xzzzz_xxxxyz[k] = -g_x_0_xzzz_xxxxyz[k] * cd_z[k] + g_x_0_xzzz_xxxxyzz[k];

                g_x_0_xzzzz_xxxxzz[k] = -g_x_0_xzzz_xxxxzz[k] * cd_z[k] + g_x_0_xzzz_xxxxzzz[k];

                g_x_0_xzzzz_xxxyyy[k] = -g_x_0_xzzz_xxxyyy[k] * cd_z[k] + g_x_0_xzzz_xxxyyyz[k];

                g_x_0_xzzzz_xxxyyz[k] = -g_x_0_xzzz_xxxyyz[k] * cd_z[k] + g_x_0_xzzz_xxxyyzz[k];

                g_x_0_xzzzz_xxxyzz[k] = -g_x_0_xzzz_xxxyzz[k] * cd_z[k] + g_x_0_xzzz_xxxyzzz[k];

                g_x_0_xzzzz_xxxzzz[k] = -g_x_0_xzzz_xxxzzz[k] * cd_z[k] + g_x_0_xzzz_xxxzzzz[k];

                g_x_0_xzzzz_xxyyyy[k] = -g_x_0_xzzz_xxyyyy[k] * cd_z[k] + g_x_0_xzzz_xxyyyyz[k];

                g_x_0_xzzzz_xxyyyz[k] = -g_x_0_xzzz_xxyyyz[k] * cd_z[k] + g_x_0_xzzz_xxyyyzz[k];

                g_x_0_xzzzz_xxyyzz[k] = -g_x_0_xzzz_xxyyzz[k] * cd_z[k] + g_x_0_xzzz_xxyyzzz[k];

                g_x_0_xzzzz_xxyzzz[k] = -g_x_0_xzzz_xxyzzz[k] * cd_z[k] + g_x_0_xzzz_xxyzzzz[k];

                g_x_0_xzzzz_xxzzzz[k] = -g_x_0_xzzz_xxzzzz[k] * cd_z[k] + g_x_0_xzzz_xxzzzzz[k];

                g_x_0_xzzzz_xyyyyy[k] = -g_x_0_xzzz_xyyyyy[k] * cd_z[k] + g_x_0_xzzz_xyyyyyz[k];

                g_x_0_xzzzz_xyyyyz[k] = -g_x_0_xzzz_xyyyyz[k] * cd_z[k] + g_x_0_xzzz_xyyyyzz[k];

                g_x_0_xzzzz_xyyyzz[k] = -g_x_0_xzzz_xyyyzz[k] * cd_z[k] + g_x_0_xzzz_xyyyzzz[k];

                g_x_0_xzzzz_xyyzzz[k] = -g_x_0_xzzz_xyyzzz[k] * cd_z[k] + g_x_0_xzzz_xyyzzzz[k];

                g_x_0_xzzzz_xyzzzz[k] = -g_x_0_xzzz_xyzzzz[k] * cd_z[k] + g_x_0_xzzz_xyzzzzz[k];

                g_x_0_xzzzz_xzzzzz[k] = -g_x_0_xzzz_xzzzzz[k] * cd_z[k] + g_x_0_xzzz_xzzzzzz[k];

                g_x_0_xzzzz_yyyyyy[k] = -g_x_0_xzzz_yyyyyy[k] * cd_z[k] + g_x_0_xzzz_yyyyyyz[k];

                g_x_0_xzzzz_yyyyyz[k] = -g_x_0_xzzz_yyyyyz[k] * cd_z[k] + g_x_0_xzzz_yyyyyzz[k];

                g_x_0_xzzzz_yyyyzz[k] = -g_x_0_xzzz_yyyyzz[k] * cd_z[k] + g_x_0_xzzz_yyyyzzz[k];

                g_x_0_xzzzz_yyyzzz[k] = -g_x_0_xzzz_yyyzzz[k] * cd_z[k] + g_x_0_xzzz_yyyzzzz[k];

                g_x_0_xzzzz_yyzzzz[k] = -g_x_0_xzzz_yyzzzz[k] * cd_z[k] + g_x_0_xzzz_yyzzzzz[k];

                g_x_0_xzzzz_yzzzzz[k] = -g_x_0_xzzz_yzzzzz[k] * cd_z[k] + g_x_0_xzzz_yzzzzzz[k];

                g_x_0_xzzzz_zzzzzz[k] = -g_x_0_xzzz_zzzzzz[k] * cd_z[k] + g_x_0_xzzz_zzzzzzz[k];
            }

            /// Set up 420-448 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 420);

            auto g_x_0_yyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 421);

            auto g_x_0_yyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 422);

            auto g_x_0_yyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 423);

            auto g_x_0_yyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 424);

            auto g_x_0_yyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 425);

            auto g_x_0_yyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 426);

            auto g_x_0_yyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 427);

            auto g_x_0_yyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 428);

            auto g_x_0_yyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 429);

            auto g_x_0_yyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 430);

            auto g_x_0_yyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 431);

            auto g_x_0_yyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 432);

            auto g_x_0_yyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 433);

            auto g_x_0_yyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 434);

            auto g_x_0_yyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 435);

            auto g_x_0_yyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 436);

            auto g_x_0_yyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 437);

            auto g_x_0_yyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 438);

            auto g_x_0_yyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 439);

            auto g_x_0_yyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 440);

            auto g_x_0_yyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 441);

            auto g_x_0_yyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 442);

            auto g_x_0_yyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 443);

            auto g_x_0_yyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 444);

            auto g_x_0_yyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 445);

            auto g_x_0_yyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 446);

            auto g_x_0_yyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 447);

            #pragma omp simd aligned(cd_y, g_x_0_yyyy_xxxxxx, g_x_0_yyyy_xxxxxxy, g_x_0_yyyy_xxxxxy, g_x_0_yyyy_xxxxxyy, g_x_0_yyyy_xxxxxyz, g_x_0_yyyy_xxxxxz, g_x_0_yyyy_xxxxyy, g_x_0_yyyy_xxxxyyy, g_x_0_yyyy_xxxxyyz, g_x_0_yyyy_xxxxyz, g_x_0_yyyy_xxxxyzz, g_x_0_yyyy_xxxxzz, g_x_0_yyyy_xxxyyy, g_x_0_yyyy_xxxyyyy, g_x_0_yyyy_xxxyyyz, g_x_0_yyyy_xxxyyz, g_x_0_yyyy_xxxyyzz, g_x_0_yyyy_xxxyzz, g_x_0_yyyy_xxxyzzz, g_x_0_yyyy_xxxzzz, g_x_0_yyyy_xxyyyy, g_x_0_yyyy_xxyyyyy, g_x_0_yyyy_xxyyyyz, g_x_0_yyyy_xxyyyz, g_x_0_yyyy_xxyyyzz, g_x_0_yyyy_xxyyzz, g_x_0_yyyy_xxyyzzz, g_x_0_yyyy_xxyzzz, g_x_0_yyyy_xxyzzzz, g_x_0_yyyy_xxzzzz, g_x_0_yyyy_xyyyyy, g_x_0_yyyy_xyyyyyy, g_x_0_yyyy_xyyyyyz, g_x_0_yyyy_xyyyyz, g_x_0_yyyy_xyyyyzz, g_x_0_yyyy_xyyyzz, g_x_0_yyyy_xyyyzzz, g_x_0_yyyy_xyyzzz, g_x_0_yyyy_xyyzzzz, g_x_0_yyyy_xyzzzz, g_x_0_yyyy_xyzzzzz, g_x_0_yyyy_xzzzzz, g_x_0_yyyy_yyyyyy, g_x_0_yyyy_yyyyyyy, g_x_0_yyyy_yyyyyyz, g_x_0_yyyy_yyyyyz, g_x_0_yyyy_yyyyyzz, g_x_0_yyyy_yyyyzz, g_x_0_yyyy_yyyyzzz, g_x_0_yyyy_yyyzzz, g_x_0_yyyy_yyyzzzz, g_x_0_yyyy_yyzzzz, g_x_0_yyyy_yyzzzzz, g_x_0_yyyy_yzzzzz, g_x_0_yyyy_yzzzzzz, g_x_0_yyyy_zzzzzz, g_x_0_yyyyy_xxxxxx, g_x_0_yyyyy_xxxxxy, g_x_0_yyyyy_xxxxxz, g_x_0_yyyyy_xxxxyy, g_x_0_yyyyy_xxxxyz, g_x_0_yyyyy_xxxxzz, g_x_0_yyyyy_xxxyyy, g_x_0_yyyyy_xxxyyz, g_x_0_yyyyy_xxxyzz, g_x_0_yyyyy_xxxzzz, g_x_0_yyyyy_xxyyyy, g_x_0_yyyyy_xxyyyz, g_x_0_yyyyy_xxyyzz, g_x_0_yyyyy_xxyzzz, g_x_0_yyyyy_xxzzzz, g_x_0_yyyyy_xyyyyy, g_x_0_yyyyy_xyyyyz, g_x_0_yyyyy_xyyyzz, g_x_0_yyyyy_xyyzzz, g_x_0_yyyyy_xyzzzz, g_x_0_yyyyy_xzzzzz, g_x_0_yyyyy_yyyyyy, g_x_0_yyyyy_yyyyyz, g_x_0_yyyyy_yyyyzz, g_x_0_yyyyy_yyyzzz, g_x_0_yyyyy_yyzzzz, g_x_0_yyyyy_yzzzzz, g_x_0_yyyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyy_xxxxxx[k] = -g_x_0_yyyy_xxxxxx[k] * cd_y[k] + g_x_0_yyyy_xxxxxxy[k];

                g_x_0_yyyyy_xxxxxy[k] = -g_x_0_yyyy_xxxxxy[k] * cd_y[k] + g_x_0_yyyy_xxxxxyy[k];

                g_x_0_yyyyy_xxxxxz[k] = -g_x_0_yyyy_xxxxxz[k] * cd_y[k] + g_x_0_yyyy_xxxxxyz[k];

                g_x_0_yyyyy_xxxxyy[k] = -g_x_0_yyyy_xxxxyy[k] * cd_y[k] + g_x_0_yyyy_xxxxyyy[k];

                g_x_0_yyyyy_xxxxyz[k] = -g_x_0_yyyy_xxxxyz[k] * cd_y[k] + g_x_0_yyyy_xxxxyyz[k];

                g_x_0_yyyyy_xxxxzz[k] = -g_x_0_yyyy_xxxxzz[k] * cd_y[k] + g_x_0_yyyy_xxxxyzz[k];

                g_x_0_yyyyy_xxxyyy[k] = -g_x_0_yyyy_xxxyyy[k] * cd_y[k] + g_x_0_yyyy_xxxyyyy[k];

                g_x_0_yyyyy_xxxyyz[k] = -g_x_0_yyyy_xxxyyz[k] * cd_y[k] + g_x_0_yyyy_xxxyyyz[k];

                g_x_0_yyyyy_xxxyzz[k] = -g_x_0_yyyy_xxxyzz[k] * cd_y[k] + g_x_0_yyyy_xxxyyzz[k];

                g_x_0_yyyyy_xxxzzz[k] = -g_x_0_yyyy_xxxzzz[k] * cd_y[k] + g_x_0_yyyy_xxxyzzz[k];

                g_x_0_yyyyy_xxyyyy[k] = -g_x_0_yyyy_xxyyyy[k] * cd_y[k] + g_x_0_yyyy_xxyyyyy[k];

                g_x_0_yyyyy_xxyyyz[k] = -g_x_0_yyyy_xxyyyz[k] * cd_y[k] + g_x_0_yyyy_xxyyyyz[k];

                g_x_0_yyyyy_xxyyzz[k] = -g_x_0_yyyy_xxyyzz[k] * cd_y[k] + g_x_0_yyyy_xxyyyzz[k];

                g_x_0_yyyyy_xxyzzz[k] = -g_x_0_yyyy_xxyzzz[k] * cd_y[k] + g_x_0_yyyy_xxyyzzz[k];

                g_x_0_yyyyy_xxzzzz[k] = -g_x_0_yyyy_xxzzzz[k] * cd_y[k] + g_x_0_yyyy_xxyzzzz[k];

                g_x_0_yyyyy_xyyyyy[k] = -g_x_0_yyyy_xyyyyy[k] * cd_y[k] + g_x_0_yyyy_xyyyyyy[k];

                g_x_0_yyyyy_xyyyyz[k] = -g_x_0_yyyy_xyyyyz[k] * cd_y[k] + g_x_0_yyyy_xyyyyyz[k];

                g_x_0_yyyyy_xyyyzz[k] = -g_x_0_yyyy_xyyyzz[k] * cd_y[k] + g_x_0_yyyy_xyyyyzz[k];

                g_x_0_yyyyy_xyyzzz[k] = -g_x_0_yyyy_xyyzzz[k] * cd_y[k] + g_x_0_yyyy_xyyyzzz[k];

                g_x_0_yyyyy_xyzzzz[k] = -g_x_0_yyyy_xyzzzz[k] * cd_y[k] + g_x_0_yyyy_xyyzzzz[k];

                g_x_0_yyyyy_xzzzzz[k] = -g_x_0_yyyy_xzzzzz[k] * cd_y[k] + g_x_0_yyyy_xyzzzzz[k];

                g_x_0_yyyyy_yyyyyy[k] = -g_x_0_yyyy_yyyyyy[k] * cd_y[k] + g_x_0_yyyy_yyyyyyy[k];

                g_x_0_yyyyy_yyyyyz[k] = -g_x_0_yyyy_yyyyyz[k] * cd_y[k] + g_x_0_yyyy_yyyyyyz[k];

                g_x_0_yyyyy_yyyyzz[k] = -g_x_0_yyyy_yyyyzz[k] * cd_y[k] + g_x_0_yyyy_yyyyyzz[k];

                g_x_0_yyyyy_yyyzzz[k] = -g_x_0_yyyy_yyyzzz[k] * cd_y[k] + g_x_0_yyyy_yyyyzzz[k];

                g_x_0_yyyyy_yyzzzz[k] = -g_x_0_yyyy_yyzzzz[k] * cd_y[k] + g_x_0_yyyy_yyyzzzz[k];

                g_x_0_yyyyy_yzzzzz[k] = -g_x_0_yyyy_yzzzzz[k] * cd_y[k] + g_x_0_yyyy_yyzzzzz[k];

                g_x_0_yyyyy_zzzzzz[k] = -g_x_0_yyyy_zzzzzz[k] * cd_y[k] + g_x_0_yyyy_yzzzzzz[k];
            }

            /// Set up 448-476 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 448);

            auto g_x_0_yyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 449);

            auto g_x_0_yyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 450);

            auto g_x_0_yyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 451);

            auto g_x_0_yyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 452);

            auto g_x_0_yyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 453);

            auto g_x_0_yyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 454);

            auto g_x_0_yyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 455);

            auto g_x_0_yyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 456);

            auto g_x_0_yyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 457);

            auto g_x_0_yyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 458);

            auto g_x_0_yyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 459);

            auto g_x_0_yyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 460);

            auto g_x_0_yyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 461);

            auto g_x_0_yyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 462);

            auto g_x_0_yyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 463);

            auto g_x_0_yyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 464);

            auto g_x_0_yyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 465);

            auto g_x_0_yyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 466);

            auto g_x_0_yyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 467);

            auto g_x_0_yyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 468);

            auto g_x_0_yyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 469);

            auto g_x_0_yyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 470);

            auto g_x_0_yyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 471);

            auto g_x_0_yyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 472);

            auto g_x_0_yyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 473);

            auto g_x_0_yyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 474);

            auto g_x_0_yyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 475);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyz_xxxxxx, g_x_0_yyyyz_xxxxxy, g_x_0_yyyyz_xxxxxz, g_x_0_yyyyz_xxxxyy, g_x_0_yyyyz_xxxxyz, g_x_0_yyyyz_xxxxzz, g_x_0_yyyyz_xxxyyy, g_x_0_yyyyz_xxxyyz, g_x_0_yyyyz_xxxyzz, g_x_0_yyyyz_xxxzzz, g_x_0_yyyyz_xxyyyy, g_x_0_yyyyz_xxyyyz, g_x_0_yyyyz_xxyyzz, g_x_0_yyyyz_xxyzzz, g_x_0_yyyyz_xxzzzz, g_x_0_yyyyz_xyyyyy, g_x_0_yyyyz_xyyyyz, g_x_0_yyyyz_xyyyzz, g_x_0_yyyyz_xyyzzz, g_x_0_yyyyz_xyzzzz, g_x_0_yyyyz_xzzzzz, g_x_0_yyyyz_yyyyyy, g_x_0_yyyyz_yyyyyz, g_x_0_yyyyz_yyyyzz, g_x_0_yyyyz_yyyzzz, g_x_0_yyyyz_yyzzzz, g_x_0_yyyyz_yzzzzz, g_x_0_yyyyz_zzzzzz, g_x_0_yyyz_xxxxxx, g_x_0_yyyz_xxxxxxy, g_x_0_yyyz_xxxxxy, g_x_0_yyyz_xxxxxyy, g_x_0_yyyz_xxxxxyz, g_x_0_yyyz_xxxxxz, g_x_0_yyyz_xxxxyy, g_x_0_yyyz_xxxxyyy, g_x_0_yyyz_xxxxyyz, g_x_0_yyyz_xxxxyz, g_x_0_yyyz_xxxxyzz, g_x_0_yyyz_xxxxzz, g_x_0_yyyz_xxxyyy, g_x_0_yyyz_xxxyyyy, g_x_0_yyyz_xxxyyyz, g_x_0_yyyz_xxxyyz, g_x_0_yyyz_xxxyyzz, g_x_0_yyyz_xxxyzz, g_x_0_yyyz_xxxyzzz, g_x_0_yyyz_xxxzzz, g_x_0_yyyz_xxyyyy, g_x_0_yyyz_xxyyyyy, g_x_0_yyyz_xxyyyyz, g_x_0_yyyz_xxyyyz, g_x_0_yyyz_xxyyyzz, g_x_0_yyyz_xxyyzz, g_x_0_yyyz_xxyyzzz, g_x_0_yyyz_xxyzzz, g_x_0_yyyz_xxyzzzz, g_x_0_yyyz_xxzzzz, g_x_0_yyyz_xyyyyy, g_x_0_yyyz_xyyyyyy, g_x_0_yyyz_xyyyyyz, g_x_0_yyyz_xyyyyz, g_x_0_yyyz_xyyyyzz, g_x_0_yyyz_xyyyzz, g_x_0_yyyz_xyyyzzz, g_x_0_yyyz_xyyzzz, g_x_0_yyyz_xyyzzzz, g_x_0_yyyz_xyzzzz, g_x_0_yyyz_xyzzzzz, g_x_0_yyyz_xzzzzz, g_x_0_yyyz_yyyyyy, g_x_0_yyyz_yyyyyyy, g_x_0_yyyz_yyyyyyz, g_x_0_yyyz_yyyyyz, g_x_0_yyyz_yyyyyzz, g_x_0_yyyz_yyyyzz, g_x_0_yyyz_yyyyzzz, g_x_0_yyyz_yyyzzz, g_x_0_yyyz_yyyzzzz, g_x_0_yyyz_yyzzzz, g_x_0_yyyz_yyzzzzz, g_x_0_yyyz_yzzzzz, g_x_0_yyyz_yzzzzzz, g_x_0_yyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyz_xxxxxx[k] = -g_x_0_yyyz_xxxxxx[k] * cd_y[k] + g_x_0_yyyz_xxxxxxy[k];

                g_x_0_yyyyz_xxxxxy[k] = -g_x_0_yyyz_xxxxxy[k] * cd_y[k] + g_x_0_yyyz_xxxxxyy[k];

                g_x_0_yyyyz_xxxxxz[k] = -g_x_0_yyyz_xxxxxz[k] * cd_y[k] + g_x_0_yyyz_xxxxxyz[k];

                g_x_0_yyyyz_xxxxyy[k] = -g_x_0_yyyz_xxxxyy[k] * cd_y[k] + g_x_0_yyyz_xxxxyyy[k];

                g_x_0_yyyyz_xxxxyz[k] = -g_x_0_yyyz_xxxxyz[k] * cd_y[k] + g_x_0_yyyz_xxxxyyz[k];

                g_x_0_yyyyz_xxxxzz[k] = -g_x_0_yyyz_xxxxzz[k] * cd_y[k] + g_x_0_yyyz_xxxxyzz[k];

                g_x_0_yyyyz_xxxyyy[k] = -g_x_0_yyyz_xxxyyy[k] * cd_y[k] + g_x_0_yyyz_xxxyyyy[k];

                g_x_0_yyyyz_xxxyyz[k] = -g_x_0_yyyz_xxxyyz[k] * cd_y[k] + g_x_0_yyyz_xxxyyyz[k];

                g_x_0_yyyyz_xxxyzz[k] = -g_x_0_yyyz_xxxyzz[k] * cd_y[k] + g_x_0_yyyz_xxxyyzz[k];

                g_x_0_yyyyz_xxxzzz[k] = -g_x_0_yyyz_xxxzzz[k] * cd_y[k] + g_x_0_yyyz_xxxyzzz[k];

                g_x_0_yyyyz_xxyyyy[k] = -g_x_0_yyyz_xxyyyy[k] * cd_y[k] + g_x_0_yyyz_xxyyyyy[k];

                g_x_0_yyyyz_xxyyyz[k] = -g_x_0_yyyz_xxyyyz[k] * cd_y[k] + g_x_0_yyyz_xxyyyyz[k];

                g_x_0_yyyyz_xxyyzz[k] = -g_x_0_yyyz_xxyyzz[k] * cd_y[k] + g_x_0_yyyz_xxyyyzz[k];

                g_x_0_yyyyz_xxyzzz[k] = -g_x_0_yyyz_xxyzzz[k] * cd_y[k] + g_x_0_yyyz_xxyyzzz[k];

                g_x_0_yyyyz_xxzzzz[k] = -g_x_0_yyyz_xxzzzz[k] * cd_y[k] + g_x_0_yyyz_xxyzzzz[k];

                g_x_0_yyyyz_xyyyyy[k] = -g_x_0_yyyz_xyyyyy[k] * cd_y[k] + g_x_0_yyyz_xyyyyyy[k];

                g_x_0_yyyyz_xyyyyz[k] = -g_x_0_yyyz_xyyyyz[k] * cd_y[k] + g_x_0_yyyz_xyyyyyz[k];

                g_x_0_yyyyz_xyyyzz[k] = -g_x_0_yyyz_xyyyzz[k] * cd_y[k] + g_x_0_yyyz_xyyyyzz[k];

                g_x_0_yyyyz_xyyzzz[k] = -g_x_0_yyyz_xyyzzz[k] * cd_y[k] + g_x_0_yyyz_xyyyzzz[k];

                g_x_0_yyyyz_xyzzzz[k] = -g_x_0_yyyz_xyzzzz[k] * cd_y[k] + g_x_0_yyyz_xyyzzzz[k];

                g_x_0_yyyyz_xzzzzz[k] = -g_x_0_yyyz_xzzzzz[k] * cd_y[k] + g_x_0_yyyz_xyzzzzz[k];

                g_x_0_yyyyz_yyyyyy[k] = -g_x_0_yyyz_yyyyyy[k] * cd_y[k] + g_x_0_yyyz_yyyyyyy[k];

                g_x_0_yyyyz_yyyyyz[k] = -g_x_0_yyyz_yyyyyz[k] * cd_y[k] + g_x_0_yyyz_yyyyyyz[k];

                g_x_0_yyyyz_yyyyzz[k] = -g_x_0_yyyz_yyyyzz[k] * cd_y[k] + g_x_0_yyyz_yyyyyzz[k];

                g_x_0_yyyyz_yyyzzz[k] = -g_x_0_yyyz_yyyzzz[k] * cd_y[k] + g_x_0_yyyz_yyyyzzz[k];

                g_x_0_yyyyz_yyzzzz[k] = -g_x_0_yyyz_yyzzzz[k] * cd_y[k] + g_x_0_yyyz_yyyzzzz[k];

                g_x_0_yyyyz_yzzzzz[k] = -g_x_0_yyyz_yzzzzz[k] * cd_y[k] + g_x_0_yyyz_yyzzzzz[k];

                g_x_0_yyyyz_zzzzzz[k] = -g_x_0_yyyz_zzzzzz[k] * cd_y[k] + g_x_0_yyyz_yzzzzzz[k];
            }

            /// Set up 476-504 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 476);

            auto g_x_0_yyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 477);

            auto g_x_0_yyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 478);

            auto g_x_0_yyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 479);

            auto g_x_0_yyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 480);

            auto g_x_0_yyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 481);

            auto g_x_0_yyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 482);

            auto g_x_0_yyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 483);

            auto g_x_0_yyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 484);

            auto g_x_0_yyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 485);

            auto g_x_0_yyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 486);

            auto g_x_0_yyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 487);

            auto g_x_0_yyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 488);

            auto g_x_0_yyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 489);

            auto g_x_0_yyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 490);

            auto g_x_0_yyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 491);

            auto g_x_0_yyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 492);

            auto g_x_0_yyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 493);

            auto g_x_0_yyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 494);

            auto g_x_0_yyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 495);

            auto g_x_0_yyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 496);

            auto g_x_0_yyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 497);

            auto g_x_0_yyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 498);

            auto g_x_0_yyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 499);

            auto g_x_0_yyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 500);

            auto g_x_0_yyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 501);

            auto g_x_0_yyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 502);

            auto g_x_0_yyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 503);

            #pragma omp simd aligned(cd_y, g_x_0_yyyzz_xxxxxx, g_x_0_yyyzz_xxxxxy, g_x_0_yyyzz_xxxxxz, g_x_0_yyyzz_xxxxyy, g_x_0_yyyzz_xxxxyz, g_x_0_yyyzz_xxxxzz, g_x_0_yyyzz_xxxyyy, g_x_0_yyyzz_xxxyyz, g_x_0_yyyzz_xxxyzz, g_x_0_yyyzz_xxxzzz, g_x_0_yyyzz_xxyyyy, g_x_0_yyyzz_xxyyyz, g_x_0_yyyzz_xxyyzz, g_x_0_yyyzz_xxyzzz, g_x_0_yyyzz_xxzzzz, g_x_0_yyyzz_xyyyyy, g_x_0_yyyzz_xyyyyz, g_x_0_yyyzz_xyyyzz, g_x_0_yyyzz_xyyzzz, g_x_0_yyyzz_xyzzzz, g_x_0_yyyzz_xzzzzz, g_x_0_yyyzz_yyyyyy, g_x_0_yyyzz_yyyyyz, g_x_0_yyyzz_yyyyzz, g_x_0_yyyzz_yyyzzz, g_x_0_yyyzz_yyzzzz, g_x_0_yyyzz_yzzzzz, g_x_0_yyyzz_zzzzzz, g_x_0_yyzz_xxxxxx, g_x_0_yyzz_xxxxxxy, g_x_0_yyzz_xxxxxy, g_x_0_yyzz_xxxxxyy, g_x_0_yyzz_xxxxxyz, g_x_0_yyzz_xxxxxz, g_x_0_yyzz_xxxxyy, g_x_0_yyzz_xxxxyyy, g_x_0_yyzz_xxxxyyz, g_x_0_yyzz_xxxxyz, g_x_0_yyzz_xxxxyzz, g_x_0_yyzz_xxxxzz, g_x_0_yyzz_xxxyyy, g_x_0_yyzz_xxxyyyy, g_x_0_yyzz_xxxyyyz, g_x_0_yyzz_xxxyyz, g_x_0_yyzz_xxxyyzz, g_x_0_yyzz_xxxyzz, g_x_0_yyzz_xxxyzzz, g_x_0_yyzz_xxxzzz, g_x_0_yyzz_xxyyyy, g_x_0_yyzz_xxyyyyy, g_x_0_yyzz_xxyyyyz, g_x_0_yyzz_xxyyyz, g_x_0_yyzz_xxyyyzz, g_x_0_yyzz_xxyyzz, g_x_0_yyzz_xxyyzzz, g_x_0_yyzz_xxyzzz, g_x_0_yyzz_xxyzzzz, g_x_0_yyzz_xxzzzz, g_x_0_yyzz_xyyyyy, g_x_0_yyzz_xyyyyyy, g_x_0_yyzz_xyyyyyz, g_x_0_yyzz_xyyyyz, g_x_0_yyzz_xyyyyzz, g_x_0_yyzz_xyyyzz, g_x_0_yyzz_xyyyzzz, g_x_0_yyzz_xyyzzz, g_x_0_yyzz_xyyzzzz, g_x_0_yyzz_xyzzzz, g_x_0_yyzz_xyzzzzz, g_x_0_yyzz_xzzzzz, g_x_0_yyzz_yyyyyy, g_x_0_yyzz_yyyyyyy, g_x_0_yyzz_yyyyyyz, g_x_0_yyzz_yyyyyz, g_x_0_yyzz_yyyyyzz, g_x_0_yyzz_yyyyzz, g_x_0_yyzz_yyyyzzz, g_x_0_yyzz_yyyzzz, g_x_0_yyzz_yyyzzzz, g_x_0_yyzz_yyzzzz, g_x_0_yyzz_yyzzzzz, g_x_0_yyzz_yzzzzz, g_x_0_yyzz_yzzzzzz, g_x_0_yyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzz_xxxxxx[k] = -g_x_0_yyzz_xxxxxx[k] * cd_y[k] + g_x_0_yyzz_xxxxxxy[k];

                g_x_0_yyyzz_xxxxxy[k] = -g_x_0_yyzz_xxxxxy[k] * cd_y[k] + g_x_0_yyzz_xxxxxyy[k];

                g_x_0_yyyzz_xxxxxz[k] = -g_x_0_yyzz_xxxxxz[k] * cd_y[k] + g_x_0_yyzz_xxxxxyz[k];

                g_x_0_yyyzz_xxxxyy[k] = -g_x_0_yyzz_xxxxyy[k] * cd_y[k] + g_x_0_yyzz_xxxxyyy[k];

                g_x_0_yyyzz_xxxxyz[k] = -g_x_0_yyzz_xxxxyz[k] * cd_y[k] + g_x_0_yyzz_xxxxyyz[k];

                g_x_0_yyyzz_xxxxzz[k] = -g_x_0_yyzz_xxxxzz[k] * cd_y[k] + g_x_0_yyzz_xxxxyzz[k];

                g_x_0_yyyzz_xxxyyy[k] = -g_x_0_yyzz_xxxyyy[k] * cd_y[k] + g_x_0_yyzz_xxxyyyy[k];

                g_x_0_yyyzz_xxxyyz[k] = -g_x_0_yyzz_xxxyyz[k] * cd_y[k] + g_x_0_yyzz_xxxyyyz[k];

                g_x_0_yyyzz_xxxyzz[k] = -g_x_0_yyzz_xxxyzz[k] * cd_y[k] + g_x_0_yyzz_xxxyyzz[k];

                g_x_0_yyyzz_xxxzzz[k] = -g_x_0_yyzz_xxxzzz[k] * cd_y[k] + g_x_0_yyzz_xxxyzzz[k];

                g_x_0_yyyzz_xxyyyy[k] = -g_x_0_yyzz_xxyyyy[k] * cd_y[k] + g_x_0_yyzz_xxyyyyy[k];

                g_x_0_yyyzz_xxyyyz[k] = -g_x_0_yyzz_xxyyyz[k] * cd_y[k] + g_x_0_yyzz_xxyyyyz[k];

                g_x_0_yyyzz_xxyyzz[k] = -g_x_0_yyzz_xxyyzz[k] * cd_y[k] + g_x_0_yyzz_xxyyyzz[k];

                g_x_0_yyyzz_xxyzzz[k] = -g_x_0_yyzz_xxyzzz[k] * cd_y[k] + g_x_0_yyzz_xxyyzzz[k];

                g_x_0_yyyzz_xxzzzz[k] = -g_x_0_yyzz_xxzzzz[k] * cd_y[k] + g_x_0_yyzz_xxyzzzz[k];

                g_x_0_yyyzz_xyyyyy[k] = -g_x_0_yyzz_xyyyyy[k] * cd_y[k] + g_x_0_yyzz_xyyyyyy[k];

                g_x_0_yyyzz_xyyyyz[k] = -g_x_0_yyzz_xyyyyz[k] * cd_y[k] + g_x_0_yyzz_xyyyyyz[k];

                g_x_0_yyyzz_xyyyzz[k] = -g_x_0_yyzz_xyyyzz[k] * cd_y[k] + g_x_0_yyzz_xyyyyzz[k];

                g_x_0_yyyzz_xyyzzz[k] = -g_x_0_yyzz_xyyzzz[k] * cd_y[k] + g_x_0_yyzz_xyyyzzz[k];

                g_x_0_yyyzz_xyzzzz[k] = -g_x_0_yyzz_xyzzzz[k] * cd_y[k] + g_x_0_yyzz_xyyzzzz[k];

                g_x_0_yyyzz_xzzzzz[k] = -g_x_0_yyzz_xzzzzz[k] * cd_y[k] + g_x_0_yyzz_xyzzzzz[k];

                g_x_0_yyyzz_yyyyyy[k] = -g_x_0_yyzz_yyyyyy[k] * cd_y[k] + g_x_0_yyzz_yyyyyyy[k];

                g_x_0_yyyzz_yyyyyz[k] = -g_x_0_yyzz_yyyyyz[k] * cd_y[k] + g_x_0_yyzz_yyyyyyz[k];

                g_x_0_yyyzz_yyyyzz[k] = -g_x_0_yyzz_yyyyzz[k] * cd_y[k] + g_x_0_yyzz_yyyyyzz[k];

                g_x_0_yyyzz_yyyzzz[k] = -g_x_0_yyzz_yyyzzz[k] * cd_y[k] + g_x_0_yyzz_yyyyzzz[k];

                g_x_0_yyyzz_yyzzzz[k] = -g_x_0_yyzz_yyzzzz[k] * cd_y[k] + g_x_0_yyzz_yyyzzzz[k];

                g_x_0_yyyzz_yzzzzz[k] = -g_x_0_yyzz_yzzzzz[k] * cd_y[k] + g_x_0_yyzz_yyzzzzz[k];

                g_x_0_yyyzz_zzzzzz[k] = -g_x_0_yyzz_zzzzzz[k] * cd_y[k] + g_x_0_yyzz_yzzzzzz[k];
            }

            /// Set up 504-532 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 504);

            auto g_x_0_yyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 505);

            auto g_x_0_yyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 506);

            auto g_x_0_yyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 507);

            auto g_x_0_yyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 508);

            auto g_x_0_yyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 509);

            auto g_x_0_yyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 510);

            auto g_x_0_yyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 511);

            auto g_x_0_yyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 512);

            auto g_x_0_yyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 513);

            auto g_x_0_yyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 514);

            auto g_x_0_yyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 515);

            auto g_x_0_yyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 516);

            auto g_x_0_yyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 517);

            auto g_x_0_yyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 518);

            auto g_x_0_yyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 519);

            auto g_x_0_yyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 520);

            auto g_x_0_yyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 521);

            auto g_x_0_yyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 522);

            auto g_x_0_yyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 523);

            auto g_x_0_yyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 524);

            auto g_x_0_yyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 525);

            auto g_x_0_yyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 526);

            auto g_x_0_yyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 527);

            auto g_x_0_yyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 528);

            auto g_x_0_yyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 529);

            auto g_x_0_yyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 530);

            auto g_x_0_yyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 531);

            #pragma omp simd aligned(cd_y, g_x_0_yyzzz_xxxxxx, g_x_0_yyzzz_xxxxxy, g_x_0_yyzzz_xxxxxz, g_x_0_yyzzz_xxxxyy, g_x_0_yyzzz_xxxxyz, g_x_0_yyzzz_xxxxzz, g_x_0_yyzzz_xxxyyy, g_x_0_yyzzz_xxxyyz, g_x_0_yyzzz_xxxyzz, g_x_0_yyzzz_xxxzzz, g_x_0_yyzzz_xxyyyy, g_x_0_yyzzz_xxyyyz, g_x_0_yyzzz_xxyyzz, g_x_0_yyzzz_xxyzzz, g_x_0_yyzzz_xxzzzz, g_x_0_yyzzz_xyyyyy, g_x_0_yyzzz_xyyyyz, g_x_0_yyzzz_xyyyzz, g_x_0_yyzzz_xyyzzz, g_x_0_yyzzz_xyzzzz, g_x_0_yyzzz_xzzzzz, g_x_0_yyzzz_yyyyyy, g_x_0_yyzzz_yyyyyz, g_x_0_yyzzz_yyyyzz, g_x_0_yyzzz_yyyzzz, g_x_0_yyzzz_yyzzzz, g_x_0_yyzzz_yzzzzz, g_x_0_yyzzz_zzzzzz, g_x_0_yzzz_xxxxxx, g_x_0_yzzz_xxxxxxy, g_x_0_yzzz_xxxxxy, g_x_0_yzzz_xxxxxyy, g_x_0_yzzz_xxxxxyz, g_x_0_yzzz_xxxxxz, g_x_0_yzzz_xxxxyy, g_x_0_yzzz_xxxxyyy, g_x_0_yzzz_xxxxyyz, g_x_0_yzzz_xxxxyz, g_x_0_yzzz_xxxxyzz, g_x_0_yzzz_xxxxzz, g_x_0_yzzz_xxxyyy, g_x_0_yzzz_xxxyyyy, g_x_0_yzzz_xxxyyyz, g_x_0_yzzz_xxxyyz, g_x_0_yzzz_xxxyyzz, g_x_0_yzzz_xxxyzz, g_x_0_yzzz_xxxyzzz, g_x_0_yzzz_xxxzzz, g_x_0_yzzz_xxyyyy, g_x_0_yzzz_xxyyyyy, g_x_0_yzzz_xxyyyyz, g_x_0_yzzz_xxyyyz, g_x_0_yzzz_xxyyyzz, g_x_0_yzzz_xxyyzz, g_x_0_yzzz_xxyyzzz, g_x_0_yzzz_xxyzzz, g_x_0_yzzz_xxyzzzz, g_x_0_yzzz_xxzzzz, g_x_0_yzzz_xyyyyy, g_x_0_yzzz_xyyyyyy, g_x_0_yzzz_xyyyyyz, g_x_0_yzzz_xyyyyz, g_x_0_yzzz_xyyyyzz, g_x_0_yzzz_xyyyzz, g_x_0_yzzz_xyyyzzz, g_x_0_yzzz_xyyzzz, g_x_0_yzzz_xyyzzzz, g_x_0_yzzz_xyzzzz, g_x_0_yzzz_xyzzzzz, g_x_0_yzzz_xzzzzz, g_x_0_yzzz_yyyyyy, g_x_0_yzzz_yyyyyyy, g_x_0_yzzz_yyyyyyz, g_x_0_yzzz_yyyyyz, g_x_0_yzzz_yyyyyzz, g_x_0_yzzz_yyyyzz, g_x_0_yzzz_yyyyzzz, g_x_0_yzzz_yyyzzz, g_x_0_yzzz_yyyzzzz, g_x_0_yzzz_yyzzzz, g_x_0_yzzz_yyzzzzz, g_x_0_yzzz_yzzzzz, g_x_0_yzzz_yzzzzzz, g_x_0_yzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzz_xxxxxx[k] = -g_x_0_yzzz_xxxxxx[k] * cd_y[k] + g_x_0_yzzz_xxxxxxy[k];

                g_x_0_yyzzz_xxxxxy[k] = -g_x_0_yzzz_xxxxxy[k] * cd_y[k] + g_x_0_yzzz_xxxxxyy[k];

                g_x_0_yyzzz_xxxxxz[k] = -g_x_0_yzzz_xxxxxz[k] * cd_y[k] + g_x_0_yzzz_xxxxxyz[k];

                g_x_0_yyzzz_xxxxyy[k] = -g_x_0_yzzz_xxxxyy[k] * cd_y[k] + g_x_0_yzzz_xxxxyyy[k];

                g_x_0_yyzzz_xxxxyz[k] = -g_x_0_yzzz_xxxxyz[k] * cd_y[k] + g_x_0_yzzz_xxxxyyz[k];

                g_x_0_yyzzz_xxxxzz[k] = -g_x_0_yzzz_xxxxzz[k] * cd_y[k] + g_x_0_yzzz_xxxxyzz[k];

                g_x_0_yyzzz_xxxyyy[k] = -g_x_0_yzzz_xxxyyy[k] * cd_y[k] + g_x_0_yzzz_xxxyyyy[k];

                g_x_0_yyzzz_xxxyyz[k] = -g_x_0_yzzz_xxxyyz[k] * cd_y[k] + g_x_0_yzzz_xxxyyyz[k];

                g_x_0_yyzzz_xxxyzz[k] = -g_x_0_yzzz_xxxyzz[k] * cd_y[k] + g_x_0_yzzz_xxxyyzz[k];

                g_x_0_yyzzz_xxxzzz[k] = -g_x_0_yzzz_xxxzzz[k] * cd_y[k] + g_x_0_yzzz_xxxyzzz[k];

                g_x_0_yyzzz_xxyyyy[k] = -g_x_0_yzzz_xxyyyy[k] * cd_y[k] + g_x_0_yzzz_xxyyyyy[k];

                g_x_0_yyzzz_xxyyyz[k] = -g_x_0_yzzz_xxyyyz[k] * cd_y[k] + g_x_0_yzzz_xxyyyyz[k];

                g_x_0_yyzzz_xxyyzz[k] = -g_x_0_yzzz_xxyyzz[k] * cd_y[k] + g_x_0_yzzz_xxyyyzz[k];

                g_x_0_yyzzz_xxyzzz[k] = -g_x_0_yzzz_xxyzzz[k] * cd_y[k] + g_x_0_yzzz_xxyyzzz[k];

                g_x_0_yyzzz_xxzzzz[k] = -g_x_0_yzzz_xxzzzz[k] * cd_y[k] + g_x_0_yzzz_xxyzzzz[k];

                g_x_0_yyzzz_xyyyyy[k] = -g_x_0_yzzz_xyyyyy[k] * cd_y[k] + g_x_0_yzzz_xyyyyyy[k];

                g_x_0_yyzzz_xyyyyz[k] = -g_x_0_yzzz_xyyyyz[k] * cd_y[k] + g_x_0_yzzz_xyyyyyz[k];

                g_x_0_yyzzz_xyyyzz[k] = -g_x_0_yzzz_xyyyzz[k] * cd_y[k] + g_x_0_yzzz_xyyyyzz[k];

                g_x_0_yyzzz_xyyzzz[k] = -g_x_0_yzzz_xyyzzz[k] * cd_y[k] + g_x_0_yzzz_xyyyzzz[k];

                g_x_0_yyzzz_xyzzzz[k] = -g_x_0_yzzz_xyzzzz[k] * cd_y[k] + g_x_0_yzzz_xyyzzzz[k];

                g_x_0_yyzzz_xzzzzz[k] = -g_x_0_yzzz_xzzzzz[k] * cd_y[k] + g_x_0_yzzz_xyzzzzz[k];

                g_x_0_yyzzz_yyyyyy[k] = -g_x_0_yzzz_yyyyyy[k] * cd_y[k] + g_x_0_yzzz_yyyyyyy[k];

                g_x_0_yyzzz_yyyyyz[k] = -g_x_0_yzzz_yyyyyz[k] * cd_y[k] + g_x_0_yzzz_yyyyyyz[k];

                g_x_0_yyzzz_yyyyzz[k] = -g_x_0_yzzz_yyyyzz[k] * cd_y[k] + g_x_0_yzzz_yyyyyzz[k];

                g_x_0_yyzzz_yyyzzz[k] = -g_x_0_yzzz_yyyzzz[k] * cd_y[k] + g_x_0_yzzz_yyyyzzz[k];

                g_x_0_yyzzz_yyzzzz[k] = -g_x_0_yzzz_yyzzzz[k] * cd_y[k] + g_x_0_yzzz_yyyzzzz[k];

                g_x_0_yyzzz_yzzzzz[k] = -g_x_0_yzzz_yzzzzz[k] * cd_y[k] + g_x_0_yzzz_yyzzzzz[k];

                g_x_0_yyzzz_zzzzzz[k] = -g_x_0_yzzz_zzzzzz[k] * cd_y[k] + g_x_0_yzzz_yzzzzzz[k];
            }

            /// Set up 532-560 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 532);

            auto g_x_0_yzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 533);

            auto g_x_0_yzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 534);

            auto g_x_0_yzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 535);

            auto g_x_0_yzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 536);

            auto g_x_0_yzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 537);

            auto g_x_0_yzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 538);

            auto g_x_0_yzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 539);

            auto g_x_0_yzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 540);

            auto g_x_0_yzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 541);

            auto g_x_0_yzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 542);

            auto g_x_0_yzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 543);

            auto g_x_0_yzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 544);

            auto g_x_0_yzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 545);

            auto g_x_0_yzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 546);

            auto g_x_0_yzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 547);

            auto g_x_0_yzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 548);

            auto g_x_0_yzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 549);

            auto g_x_0_yzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 550);

            auto g_x_0_yzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 551);

            auto g_x_0_yzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 552);

            auto g_x_0_yzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 553);

            auto g_x_0_yzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 554);

            auto g_x_0_yzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 555);

            auto g_x_0_yzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 556);

            auto g_x_0_yzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 557);

            auto g_x_0_yzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 558);

            auto g_x_0_yzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 559);

            #pragma omp simd aligned(cd_y, g_x_0_yzzzz_xxxxxx, g_x_0_yzzzz_xxxxxy, g_x_0_yzzzz_xxxxxz, g_x_0_yzzzz_xxxxyy, g_x_0_yzzzz_xxxxyz, g_x_0_yzzzz_xxxxzz, g_x_0_yzzzz_xxxyyy, g_x_0_yzzzz_xxxyyz, g_x_0_yzzzz_xxxyzz, g_x_0_yzzzz_xxxzzz, g_x_0_yzzzz_xxyyyy, g_x_0_yzzzz_xxyyyz, g_x_0_yzzzz_xxyyzz, g_x_0_yzzzz_xxyzzz, g_x_0_yzzzz_xxzzzz, g_x_0_yzzzz_xyyyyy, g_x_0_yzzzz_xyyyyz, g_x_0_yzzzz_xyyyzz, g_x_0_yzzzz_xyyzzz, g_x_0_yzzzz_xyzzzz, g_x_0_yzzzz_xzzzzz, g_x_0_yzzzz_yyyyyy, g_x_0_yzzzz_yyyyyz, g_x_0_yzzzz_yyyyzz, g_x_0_yzzzz_yyyzzz, g_x_0_yzzzz_yyzzzz, g_x_0_yzzzz_yzzzzz, g_x_0_yzzzz_zzzzzz, g_x_0_zzzz_xxxxxx, g_x_0_zzzz_xxxxxxy, g_x_0_zzzz_xxxxxy, g_x_0_zzzz_xxxxxyy, g_x_0_zzzz_xxxxxyz, g_x_0_zzzz_xxxxxz, g_x_0_zzzz_xxxxyy, g_x_0_zzzz_xxxxyyy, g_x_0_zzzz_xxxxyyz, g_x_0_zzzz_xxxxyz, g_x_0_zzzz_xxxxyzz, g_x_0_zzzz_xxxxzz, g_x_0_zzzz_xxxyyy, g_x_0_zzzz_xxxyyyy, g_x_0_zzzz_xxxyyyz, g_x_0_zzzz_xxxyyz, g_x_0_zzzz_xxxyyzz, g_x_0_zzzz_xxxyzz, g_x_0_zzzz_xxxyzzz, g_x_0_zzzz_xxxzzz, g_x_0_zzzz_xxyyyy, g_x_0_zzzz_xxyyyyy, g_x_0_zzzz_xxyyyyz, g_x_0_zzzz_xxyyyz, g_x_0_zzzz_xxyyyzz, g_x_0_zzzz_xxyyzz, g_x_0_zzzz_xxyyzzz, g_x_0_zzzz_xxyzzz, g_x_0_zzzz_xxyzzzz, g_x_0_zzzz_xxzzzz, g_x_0_zzzz_xyyyyy, g_x_0_zzzz_xyyyyyy, g_x_0_zzzz_xyyyyyz, g_x_0_zzzz_xyyyyz, g_x_0_zzzz_xyyyyzz, g_x_0_zzzz_xyyyzz, g_x_0_zzzz_xyyyzzz, g_x_0_zzzz_xyyzzz, g_x_0_zzzz_xyyzzzz, g_x_0_zzzz_xyzzzz, g_x_0_zzzz_xyzzzzz, g_x_0_zzzz_xzzzzz, g_x_0_zzzz_yyyyyy, g_x_0_zzzz_yyyyyyy, g_x_0_zzzz_yyyyyyz, g_x_0_zzzz_yyyyyz, g_x_0_zzzz_yyyyyzz, g_x_0_zzzz_yyyyzz, g_x_0_zzzz_yyyyzzz, g_x_0_zzzz_yyyzzz, g_x_0_zzzz_yyyzzzz, g_x_0_zzzz_yyzzzz, g_x_0_zzzz_yyzzzzz, g_x_0_zzzz_yzzzzz, g_x_0_zzzz_yzzzzzz, g_x_0_zzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzz_xxxxxx[k] = -g_x_0_zzzz_xxxxxx[k] * cd_y[k] + g_x_0_zzzz_xxxxxxy[k];

                g_x_0_yzzzz_xxxxxy[k] = -g_x_0_zzzz_xxxxxy[k] * cd_y[k] + g_x_0_zzzz_xxxxxyy[k];

                g_x_0_yzzzz_xxxxxz[k] = -g_x_0_zzzz_xxxxxz[k] * cd_y[k] + g_x_0_zzzz_xxxxxyz[k];

                g_x_0_yzzzz_xxxxyy[k] = -g_x_0_zzzz_xxxxyy[k] * cd_y[k] + g_x_0_zzzz_xxxxyyy[k];

                g_x_0_yzzzz_xxxxyz[k] = -g_x_0_zzzz_xxxxyz[k] * cd_y[k] + g_x_0_zzzz_xxxxyyz[k];

                g_x_0_yzzzz_xxxxzz[k] = -g_x_0_zzzz_xxxxzz[k] * cd_y[k] + g_x_0_zzzz_xxxxyzz[k];

                g_x_0_yzzzz_xxxyyy[k] = -g_x_0_zzzz_xxxyyy[k] * cd_y[k] + g_x_0_zzzz_xxxyyyy[k];

                g_x_0_yzzzz_xxxyyz[k] = -g_x_0_zzzz_xxxyyz[k] * cd_y[k] + g_x_0_zzzz_xxxyyyz[k];

                g_x_0_yzzzz_xxxyzz[k] = -g_x_0_zzzz_xxxyzz[k] * cd_y[k] + g_x_0_zzzz_xxxyyzz[k];

                g_x_0_yzzzz_xxxzzz[k] = -g_x_0_zzzz_xxxzzz[k] * cd_y[k] + g_x_0_zzzz_xxxyzzz[k];

                g_x_0_yzzzz_xxyyyy[k] = -g_x_0_zzzz_xxyyyy[k] * cd_y[k] + g_x_0_zzzz_xxyyyyy[k];

                g_x_0_yzzzz_xxyyyz[k] = -g_x_0_zzzz_xxyyyz[k] * cd_y[k] + g_x_0_zzzz_xxyyyyz[k];

                g_x_0_yzzzz_xxyyzz[k] = -g_x_0_zzzz_xxyyzz[k] * cd_y[k] + g_x_0_zzzz_xxyyyzz[k];

                g_x_0_yzzzz_xxyzzz[k] = -g_x_0_zzzz_xxyzzz[k] * cd_y[k] + g_x_0_zzzz_xxyyzzz[k];

                g_x_0_yzzzz_xxzzzz[k] = -g_x_0_zzzz_xxzzzz[k] * cd_y[k] + g_x_0_zzzz_xxyzzzz[k];

                g_x_0_yzzzz_xyyyyy[k] = -g_x_0_zzzz_xyyyyy[k] * cd_y[k] + g_x_0_zzzz_xyyyyyy[k];

                g_x_0_yzzzz_xyyyyz[k] = -g_x_0_zzzz_xyyyyz[k] * cd_y[k] + g_x_0_zzzz_xyyyyyz[k];

                g_x_0_yzzzz_xyyyzz[k] = -g_x_0_zzzz_xyyyzz[k] * cd_y[k] + g_x_0_zzzz_xyyyyzz[k];

                g_x_0_yzzzz_xyyzzz[k] = -g_x_0_zzzz_xyyzzz[k] * cd_y[k] + g_x_0_zzzz_xyyyzzz[k];

                g_x_0_yzzzz_xyzzzz[k] = -g_x_0_zzzz_xyzzzz[k] * cd_y[k] + g_x_0_zzzz_xyyzzzz[k];

                g_x_0_yzzzz_xzzzzz[k] = -g_x_0_zzzz_xzzzzz[k] * cd_y[k] + g_x_0_zzzz_xyzzzzz[k];

                g_x_0_yzzzz_yyyyyy[k] = -g_x_0_zzzz_yyyyyy[k] * cd_y[k] + g_x_0_zzzz_yyyyyyy[k];

                g_x_0_yzzzz_yyyyyz[k] = -g_x_0_zzzz_yyyyyz[k] * cd_y[k] + g_x_0_zzzz_yyyyyyz[k];

                g_x_0_yzzzz_yyyyzz[k] = -g_x_0_zzzz_yyyyzz[k] * cd_y[k] + g_x_0_zzzz_yyyyyzz[k];

                g_x_0_yzzzz_yyyzzz[k] = -g_x_0_zzzz_yyyzzz[k] * cd_y[k] + g_x_0_zzzz_yyyyzzz[k];

                g_x_0_yzzzz_yyzzzz[k] = -g_x_0_zzzz_yyzzzz[k] * cd_y[k] + g_x_0_zzzz_yyyzzzz[k];

                g_x_0_yzzzz_yzzzzz[k] = -g_x_0_zzzz_yzzzzz[k] * cd_y[k] + g_x_0_zzzz_yyzzzzz[k];

                g_x_0_yzzzz_zzzzzz[k] = -g_x_0_zzzz_zzzzzz[k] * cd_y[k] + g_x_0_zzzz_yzzzzzz[k];
            }

            /// Set up 560-588 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 560);

            auto g_x_0_zzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 561);

            auto g_x_0_zzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 562);

            auto g_x_0_zzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 563);

            auto g_x_0_zzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 564);

            auto g_x_0_zzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 565);

            auto g_x_0_zzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 566);

            auto g_x_0_zzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 567);

            auto g_x_0_zzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 568);

            auto g_x_0_zzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 569);

            auto g_x_0_zzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 570);

            auto g_x_0_zzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 571);

            auto g_x_0_zzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 572);

            auto g_x_0_zzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 573);

            auto g_x_0_zzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 574);

            auto g_x_0_zzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 575);

            auto g_x_0_zzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 576);

            auto g_x_0_zzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 577);

            auto g_x_0_zzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 578);

            auto g_x_0_zzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 579);

            auto g_x_0_zzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 580);

            auto g_x_0_zzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 581);

            auto g_x_0_zzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 582);

            auto g_x_0_zzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 583);

            auto g_x_0_zzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 584);

            auto g_x_0_zzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 585);

            auto g_x_0_zzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 586);

            auto g_x_0_zzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 587);

            #pragma omp simd aligned(cd_z, g_x_0_zzzz_xxxxxx, g_x_0_zzzz_xxxxxxz, g_x_0_zzzz_xxxxxy, g_x_0_zzzz_xxxxxyz, g_x_0_zzzz_xxxxxz, g_x_0_zzzz_xxxxxzz, g_x_0_zzzz_xxxxyy, g_x_0_zzzz_xxxxyyz, g_x_0_zzzz_xxxxyz, g_x_0_zzzz_xxxxyzz, g_x_0_zzzz_xxxxzz, g_x_0_zzzz_xxxxzzz, g_x_0_zzzz_xxxyyy, g_x_0_zzzz_xxxyyyz, g_x_0_zzzz_xxxyyz, g_x_0_zzzz_xxxyyzz, g_x_0_zzzz_xxxyzz, g_x_0_zzzz_xxxyzzz, g_x_0_zzzz_xxxzzz, g_x_0_zzzz_xxxzzzz, g_x_0_zzzz_xxyyyy, g_x_0_zzzz_xxyyyyz, g_x_0_zzzz_xxyyyz, g_x_0_zzzz_xxyyyzz, g_x_0_zzzz_xxyyzz, g_x_0_zzzz_xxyyzzz, g_x_0_zzzz_xxyzzz, g_x_0_zzzz_xxyzzzz, g_x_0_zzzz_xxzzzz, g_x_0_zzzz_xxzzzzz, g_x_0_zzzz_xyyyyy, g_x_0_zzzz_xyyyyyz, g_x_0_zzzz_xyyyyz, g_x_0_zzzz_xyyyyzz, g_x_0_zzzz_xyyyzz, g_x_0_zzzz_xyyyzzz, g_x_0_zzzz_xyyzzz, g_x_0_zzzz_xyyzzzz, g_x_0_zzzz_xyzzzz, g_x_0_zzzz_xyzzzzz, g_x_0_zzzz_xzzzzz, g_x_0_zzzz_xzzzzzz, g_x_0_zzzz_yyyyyy, g_x_0_zzzz_yyyyyyz, g_x_0_zzzz_yyyyyz, g_x_0_zzzz_yyyyyzz, g_x_0_zzzz_yyyyzz, g_x_0_zzzz_yyyyzzz, g_x_0_zzzz_yyyzzz, g_x_0_zzzz_yyyzzzz, g_x_0_zzzz_yyzzzz, g_x_0_zzzz_yyzzzzz, g_x_0_zzzz_yzzzzz, g_x_0_zzzz_yzzzzzz, g_x_0_zzzz_zzzzzz, g_x_0_zzzz_zzzzzzz, g_x_0_zzzzz_xxxxxx, g_x_0_zzzzz_xxxxxy, g_x_0_zzzzz_xxxxxz, g_x_0_zzzzz_xxxxyy, g_x_0_zzzzz_xxxxyz, g_x_0_zzzzz_xxxxzz, g_x_0_zzzzz_xxxyyy, g_x_0_zzzzz_xxxyyz, g_x_0_zzzzz_xxxyzz, g_x_0_zzzzz_xxxzzz, g_x_0_zzzzz_xxyyyy, g_x_0_zzzzz_xxyyyz, g_x_0_zzzzz_xxyyzz, g_x_0_zzzzz_xxyzzz, g_x_0_zzzzz_xxzzzz, g_x_0_zzzzz_xyyyyy, g_x_0_zzzzz_xyyyyz, g_x_0_zzzzz_xyyyzz, g_x_0_zzzzz_xyyzzz, g_x_0_zzzzz_xyzzzz, g_x_0_zzzzz_xzzzzz, g_x_0_zzzzz_yyyyyy, g_x_0_zzzzz_yyyyyz, g_x_0_zzzzz_yyyyzz, g_x_0_zzzzz_yyyzzz, g_x_0_zzzzz_yyzzzz, g_x_0_zzzzz_yzzzzz, g_x_0_zzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzz_xxxxxx[k] = -g_x_0_zzzz_xxxxxx[k] * cd_z[k] + g_x_0_zzzz_xxxxxxz[k];

                g_x_0_zzzzz_xxxxxy[k] = -g_x_0_zzzz_xxxxxy[k] * cd_z[k] + g_x_0_zzzz_xxxxxyz[k];

                g_x_0_zzzzz_xxxxxz[k] = -g_x_0_zzzz_xxxxxz[k] * cd_z[k] + g_x_0_zzzz_xxxxxzz[k];

                g_x_0_zzzzz_xxxxyy[k] = -g_x_0_zzzz_xxxxyy[k] * cd_z[k] + g_x_0_zzzz_xxxxyyz[k];

                g_x_0_zzzzz_xxxxyz[k] = -g_x_0_zzzz_xxxxyz[k] * cd_z[k] + g_x_0_zzzz_xxxxyzz[k];

                g_x_0_zzzzz_xxxxzz[k] = -g_x_0_zzzz_xxxxzz[k] * cd_z[k] + g_x_0_zzzz_xxxxzzz[k];

                g_x_0_zzzzz_xxxyyy[k] = -g_x_0_zzzz_xxxyyy[k] * cd_z[k] + g_x_0_zzzz_xxxyyyz[k];

                g_x_0_zzzzz_xxxyyz[k] = -g_x_0_zzzz_xxxyyz[k] * cd_z[k] + g_x_0_zzzz_xxxyyzz[k];

                g_x_0_zzzzz_xxxyzz[k] = -g_x_0_zzzz_xxxyzz[k] * cd_z[k] + g_x_0_zzzz_xxxyzzz[k];

                g_x_0_zzzzz_xxxzzz[k] = -g_x_0_zzzz_xxxzzz[k] * cd_z[k] + g_x_0_zzzz_xxxzzzz[k];

                g_x_0_zzzzz_xxyyyy[k] = -g_x_0_zzzz_xxyyyy[k] * cd_z[k] + g_x_0_zzzz_xxyyyyz[k];

                g_x_0_zzzzz_xxyyyz[k] = -g_x_0_zzzz_xxyyyz[k] * cd_z[k] + g_x_0_zzzz_xxyyyzz[k];

                g_x_0_zzzzz_xxyyzz[k] = -g_x_0_zzzz_xxyyzz[k] * cd_z[k] + g_x_0_zzzz_xxyyzzz[k];

                g_x_0_zzzzz_xxyzzz[k] = -g_x_0_zzzz_xxyzzz[k] * cd_z[k] + g_x_0_zzzz_xxyzzzz[k];

                g_x_0_zzzzz_xxzzzz[k] = -g_x_0_zzzz_xxzzzz[k] * cd_z[k] + g_x_0_zzzz_xxzzzzz[k];

                g_x_0_zzzzz_xyyyyy[k] = -g_x_0_zzzz_xyyyyy[k] * cd_z[k] + g_x_0_zzzz_xyyyyyz[k];

                g_x_0_zzzzz_xyyyyz[k] = -g_x_0_zzzz_xyyyyz[k] * cd_z[k] + g_x_0_zzzz_xyyyyzz[k];

                g_x_0_zzzzz_xyyyzz[k] = -g_x_0_zzzz_xyyyzz[k] * cd_z[k] + g_x_0_zzzz_xyyyzzz[k];

                g_x_0_zzzzz_xyyzzz[k] = -g_x_0_zzzz_xyyzzz[k] * cd_z[k] + g_x_0_zzzz_xyyzzzz[k];

                g_x_0_zzzzz_xyzzzz[k] = -g_x_0_zzzz_xyzzzz[k] * cd_z[k] + g_x_0_zzzz_xyzzzzz[k];

                g_x_0_zzzzz_xzzzzz[k] = -g_x_0_zzzz_xzzzzz[k] * cd_z[k] + g_x_0_zzzz_xzzzzzz[k];

                g_x_0_zzzzz_yyyyyy[k] = -g_x_0_zzzz_yyyyyy[k] * cd_z[k] + g_x_0_zzzz_yyyyyyz[k];

                g_x_0_zzzzz_yyyyyz[k] = -g_x_0_zzzz_yyyyyz[k] * cd_z[k] + g_x_0_zzzz_yyyyyzz[k];

                g_x_0_zzzzz_yyyyzz[k] = -g_x_0_zzzz_yyyyzz[k] * cd_z[k] + g_x_0_zzzz_yyyyzzz[k];

                g_x_0_zzzzz_yyyzzz[k] = -g_x_0_zzzz_yyyzzz[k] * cd_z[k] + g_x_0_zzzz_yyyzzzz[k];

                g_x_0_zzzzz_yyzzzz[k] = -g_x_0_zzzz_yyzzzz[k] * cd_z[k] + g_x_0_zzzz_yyzzzzz[k];

                g_x_0_zzzzz_yzzzzz[k] = -g_x_0_zzzz_yzzzzz[k] * cd_z[k] + g_x_0_zzzz_yzzzzzz[k];

                g_x_0_zzzzz_zzzzzz[k] = -g_x_0_zzzz_zzzzzz[k] * cd_z[k] + g_x_0_zzzz_zzzzzzz[k];
            }
            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxx_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 0);

            auto g_y_0_xxxxx_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 1);

            auto g_y_0_xxxxx_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 2);

            auto g_y_0_xxxxx_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 3);

            auto g_y_0_xxxxx_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 4);

            auto g_y_0_xxxxx_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 5);

            auto g_y_0_xxxxx_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 6);

            auto g_y_0_xxxxx_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 7);

            auto g_y_0_xxxxx_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 8);

            auto g_y_0_xxxxx_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 9);

            auto g_y_0_xxxxx_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 10);

            auto g_y_0_xxxxx_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 11);

            auto g_y_0_xxxxx_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 12);

            auto g_y_0_xxxxx_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 13);

            auto g_y_0_xxxxx_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 14);

            auto g_y_0_xxxxx_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 15);

            auto g_y_0_xxxxx_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 16);

            auto g_y_0_xxxxx_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 17);

            auto g_y_0_xxxxx_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 18);

            auto g_y_0_xxxxx_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 19);

            auto g_y_0_xxxxx_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 20);

            auto g_y_0_xxxxx_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 21);

            auto g_y_0_xxxxx_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 22);

            auto g_y_0_xxxxx_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 23);

            auto g_y_0_xxxxx_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 24);

            auto g_y_0_xxxxx_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 25);

            auto g_y_0_xxxxx_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 26);

            auto g_y_0_xxxxx_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 27);

            #pragma omp simd aligned(cd_x, g_y_0_xxxx_xxxxxx, g_y_0_xxxx_xxxxxxx, g_y_0_xxxx_xxxxxxy, g_y_0_xxxx_xxxxxxz, g_y_0_xxxx_xxxxxy, g_y_0_xxxx_xxxxxyy, g_y_0_xxxx_xxxxxyz, g_y_0_xxxx_xxxxxz, g_y_0_xxxx_xxxxxzz, g_y_0_xxxx_xxxxyy, g_y_0_xxxx_xxxxyyy, g_y_0_xxxx_xxxxyyz, g_y_0_xxxx_xxxxyz, g_y_0_xxxx_xxxxyzz, g_y_0_xxxx_xxxxzz, g_y_0_xxxx_xxxxzzz, g_y_0_xxxx_xxxyyy, g_y_0_xxxx_xxxyyyy, g_y_0_xxxx_xxxyyyz, g_y_0_xxxx_xxxyyz, g_y_0_xxxx_xxxyyzz, g_y_0_xxxx_xxxyzz, g_y_0_xxxx_xxxyzzz, g_y_0_xxxx_xxxzzz, g_y_0_xxxx_xxxzzzz, g_y_0_xxxx_xxyyyy, g_y_0_xxxx_xxyyyyy, g_y_0_xxxx_xxyyyyz, g_y_0_xxxx_xxyyyz, g_y_0_xxxx_xxyyyzz, g_y_0_xxxx_xxyyzz, g_y_0_xxxx_xxyyzzz, g_y_0_xxxx_xxyzzz, g_y_0_xxxx_xxyzzzz, g_y_0_xxxx_xxzzzz, g_y_0_xxxx_xxzzzzz, g_y_0_xxxx_xyyyyy, g_y_0_xxxx_xyyyyyy, g_y_0_xxxx_xyyyyyz, g_y_0_xxxx_xyyyyz, g_y_0_xxxx_xyyyyzz, g_y_0_xxxx_xyyyzz, g_y_0_xxxx_xyyyzzz, g_y_0_xxxx_xyyzzz, g_y_0_xxxx_xyyzzzz, g_y_0_xxxx_xyzzzz, g_y_0_xxxx_xyzzzzz, g_y_0_xxxx_xzzzzz, g_y_0_xxxx_xzzzzzz, g_y_0_xxxx_yyyyyy, g_y_0_xxxx_yyyyyz, g_y_0_xxxx_yyyyzz, g_y_0_xxxx_yyyzzz, g_y_0_xxxx_yyzzzz, g_y_0_xxxx_yzzzzz, g_y_0_xxxx_zzzzzz, g_y_0_xxxxx_xxxxxx, g_y_0_xxxxx_xxxxxy, g_y_0_xxxxx_xxxxxz, g_y_0_xxxxx_xxxxyy, g_y_0_xxxxx_xxxxyz, g_y_0_xxxxx_xxxxzz, g_y_0_xxxxx_xxxyyy, g_y_0_xxxxx_xxxyyz, g_y_0_xxxxx_xxxyzz, g_y_0_xxxxx_xxxzzz, g_y_0_xxxxx_xxyyyy, g_y_0_xxxxx_xxyyyz, g_y_0_xxxxx_xxyyzz, g_y_0_xxxxx_xxyzzz, g_y_0_xxxxx_xxzzzz, g_y_0_xxxxx_xyyyyy, g_y_0_xxxxx_xyyyyz, g_y_0_xxxxx_xyyyzz, g_y_0_xxxxx_xyyzzz, g_y_0_xxxxx_xyzzzz, g_y_0_xxxxx_xzzzzz, g_y_0_xxxxx_yyyyyy, g_y_0_xxxxx_yyyyyz, g_y_0_xxxxx_yyyyzz, g_y_0_xxxxx_yyyzzz, g_y_0_xxxxx_yyzzzz, g_y_0_xxxxx_yzzzzz, g_y_0_xxxxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxx_xxxxxx[k] = -g_y_0_xxxx_xxxxxx[k] * cd_x[k] + g_y_0_xxxx_xxxxxxx[k];

                g_y_0_xxxxx_xxxxxy[k] = -g_y_0_xxxx_xxxxxy[k] * cd_x[k] + g_y_0_xxxx_xxxxxxy[k];

                g_y_0_xxxxx_xxxxxz[k] = -g_y_0_xxxx_xxxxxz[k] * cd_x[k] + g_y_0_xxxx_xxxxxxz[k];

                g_y_0_xxxxx_xxxxyy[k] = -g_y_0_xxxx_xxxxyy[k] * cd_x[k] + g_y_0_xxxx_xxxxxyy[k];

                g_y_0_xxxxx_xxxxyz[k] = -g_y_0_xxxx_xxxxyz[k] * cd_x[k] + g_y_0_xxxx_xxxxxyz[k];

                g_y_0_xxxxx_xxxxzz[k] = -g_y_0_xxxx_xxxxzz[k] * cd_x[k] + g_y_0_xxxx_xxxxxzz[k];

                g_y_0_xxxxx_xxxyyy[k] = -g_y_0_xxxx_xxxyyy[k] * cd_x[k] + g_y_0_xxxx_xxxxyyy[k];

                g_y_0_xxxxx_xxxyyz[k] = -g_y_0_xxxx_xxxyyz[k] * cd_x[k] + g_y_0_xxxx_xxxxyyz[k];

                g_y_0_xxxxx_xxxyzz[k] = -g_y_0_xxxx_xxxyzz[k] * cd_x[k] + g_y_0_xxxx_xxxxyzz[k];

                g_y_0_xxxxx_xxxzzz[k] = -g_y_0_xxxx_xxxzzz[k] * cd_x[k] + g_y_0_xxxx_xxxxzzz[k];

                g_y_0_xxxxx_xxyyyy[k] = -g_y_0_xxxx_xxyyyy[k] * cd_x[k] + g_y_0_xxxx_xxxyyyy[k];

                g_y_0_xxxxx_xxyyyz[k] = -g_y_0_xxxx_xxyyyz[k] * cd_x[k] + g_y_0_xxxx_xxxyyyz[k];

                g_y_0_xxxxx_xxyyzz[k] = -g_y_0_xxxx_xxyyzz[k] * cd_x[k] + g_y_0_xxxx_xxxyyzz[k];

                g_y_0_xxxxx_xxyzzz[k] = -g_y_0_xxxx_xxyzzz[k] * cd_x[k] + g_y_0_xxxx_xxxyzzz[k];

                g_y_0_xxxxx_xxzzzz[k] = -g_y_0_xxxx_xxzzzz[k] * cd_x[k] + g_y_0_xxxx_xxxzzzz[k];

                g_y_0_xxxxx_xyyyyy[k] = -g_y_0_xxxx_xyyyyy[k] * cd_x[k] + g_y_0_xxxx_xxyyyyy[k];

                g_y_0_xxxxx_xyyyyz[k] = -g_y_0_xxxx_xyyyyz[k] * cd_x[k] + g_y_0_xxxx_xxyyyyz[k];

                g_y_0_xxxxx_xyyyzz[k] = -g_y_0_xxxx_xyyyzz[k] * cd_x[k] + g_y_0_xxxx_xxyyyzz[k];

                g_y_0_xxxxx_xyyzzz[k] = -g_y_0_xxxx_xyyzzz[k] * cd_x[k] + g_y_0_xxxx_xxyyzzz[k];

                g_y_0_xxxxx_xyzzzz[k] = -g_y_0_xxxx_xyzzzz[k] * cd_x[k] + g_y_0_xxxx_xxyzzzz[k];

                g_y_0_xxxxx_xzzzzz[k] = -g_y_0_xxxx_xzzzzz[k] * cd_x[k] + g_y_0_xxxx_xxzzzzz[k];

                g_y_0_xxxxx_yyyyyy[k] = -g_y_0_xxxx_yyyyyy[k] * cd_x[k] + g_y_0_xxxx_xyyyyyy[k];

                g_y_0_xxxxx_yyyyyz[k] = -g_y_0_xxxx_yyyyyz[k] * cd_x[k] + g_y_0_xxxx_xyyyyyz[k];

                g_y_0_xxxxx_yyyyzz[k] = -g_y_0_xxxx_yyyyzz[k] * cd_x[k] + g_y_0_xxxx_xyyyyzz[k];

                g_y_0_xxxxx_yyyzzz[k] = -g_y_0_xxxx_yyyzzz[k] * cd_x[k] + g_y_0_xxxx_xyyyzzz[k];

                g_y_0_xxxxx_yyzzzz[k] = -g_y_0_xxxx_yyzzzz[k] * cd_x[k] + g_y_0_xxxx_xyyzzzz[k];

                g_y_0_xxxxx_yzzzzz[k] = -g_y_0_xxxx_yzzzzz[k] * cd_x[k] + g_y_0_xxxx_xyzzzzz[k];

                g_y_0_xxxxx_zzzzzz[k] = -g_y_0_xxxx_zzzzzz[k] * cd_x[k] + g_y_0_xxxx_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxy_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 28);

            auto g_y_0_xxxxy_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 29);

            auto g_y_0_xxxxy_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 30);

            auto g_y_0_xxxxy_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 31);

            auto g_y_0_xxxxy_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 32);

            auto g_y_0_xxxxy_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 33);

            auto g_y_0_xxxxy_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 34);

            auto g_y_0_xxxxy_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 35);

            auto g_y_0_xxxxy_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 36);

            auto g_y_0_xxxxy_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 37);

            auto g_y_0_xxxxy_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 38);

            auto g_y_0_xxxxy_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 39);

            auto g_y_0_xxxxy_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 40);

            auto g_y_0_xxxxy_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 41);

            auto g_y_0_xxxxy_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 42);

            auto g_y_0_xxxxy_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 43);

            auto g_y_0_xxxxy_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 44);

            auto g_y_0_xxxxy_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 45);

            auto g_y_0_xxxxy_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 46);

            auto g_y_0_xxxxy_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 47);

            auto g_y_0_xxxxy_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 48);

            auto g_y_0_xxxxy_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 49);

            auto g_y_0_xxxxy_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 50);

            auto g_y_0_xxxxy_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 51);

            auto g_y_0_xxxxy_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 52);

            auto g_y_0_xxxxy_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 53);

            auto g_y_0_xxxxy_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 54);

            auto g_y_0_xxxxy_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 55);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxy_xxxxxx, g_y_0_xxxxy_xxxxxy, g_y_0_xxxxy_xxxxxz, g_y_0_xxxxy_xxxxyy, g_y_0_xxxxy_xxxxyz, g_y_0_xxxxy_xxxxzz, g_y_0_xxxxy_xxxyyy, g_y_0_xxxxy_xxxyyz, g_y_0_xxxxy_xxxyzz, g_y_0_xxxxy_xxxzzz, g_y_0_xxxxy_xxyyyy, g_y_0_xxxxy_xxyyyz, g_y_0_xxxxy_xxyyzz, g_y_0_xxxxy_xxyzzz, g_y_0_xxxxy_xxzzzz, g_y_0_xxxxy_xyyyyy, g_y_0_xxxxy_xyyyyz, g_y_0_xxxxy_xyyyzz, g_y_0_xxxxy_xyyzzz, g_y_0_xxxxy_xyzzzz, g_y_0_xxxxy_xzzzzz, g_y_0_xxxxy_yyyyyy, g_y_0_xxxxy_yyyyyz, g_y_0_xxxxy_yyyyzz, g_y_0_xxxxy_yyyzzz, g_y_0_xxxxy_yyzzzz, g_y_0_xxxxy_yzzzzz, g_y_0_xxxxy_zzzzzz, g_y_0_xxxy_xxxxxx, g_y_0_xxxy_xxxxxxx, g_y_0_xxxy_xxxxxxy, g_y_0_xxxy_xxxxxxz, g_y_0_xxxy_xxxxxy, g_y_0_xxxy_xxxxxyy, g_y_0_xxxy_xxxxxyz, g_y_0_xxxy_xxxxxz, g_y_0_xxxy_xxxxxzz, g_y_0_xxxy_xxxxyy, g_y_0_xxxy_xxxxyyy, g_y_0_xxxy_xxxxyyz, g_y_0_xxxy_xxxxyz, g_y_0_xxxy_xxxxyzz, g_y_0_xxxy_xxxxzz, g_y_0_xxxy_xxxxzzz, g_y_0_xxxy_xxxyyy, g_y_0_xxxy_xxxyyyy, g_y_0_xxxy_xxxyyyz, g_y_0_xxxy_xxxyyz, g_y_0_xxxy_xxxyyzz, g_y_0_xxxy_xxxyzz, g_y_0_xxxy_xxxyzzz, g_y_0_xxxy_xxxzzz, g_y_0_xxxy_xxxzzzz, g_y_0_xxxy_xxyyyy, g_y_0_xxxy_xxyyyyy, g_y_0_xxxy_xxyyyyz, g_y_0_xxxy_xxyyyz, g_y_0_xxxy_xxyyyzz, g_y_0_xxxy_xxyyzz, g_y_0_xxxy_xxyyzzz, g_y_0_xxxy_xxyzzz, g_y_0_xxxy_xxyzzzz, g_y_0_xxxy_xxzzzz, g_y_0_xxxy_xxzzzzz, g_y_0_xxxy_xyyyyy, g_y_0_xxxy_xyyyyyy, g_y_0_xxxy_xyyyyyz, g_y_0_xxxy_xyyyyz, g_y_0_xxxy_xyyyyzz, g_y_0_xxxy_xyyyzz, g_y_0_xxxy_xyyyzzz, g_y_0_xxxy_xyyzzz, g_y_0_xxxy_xyyzzzz, g_y_0_xxxy_xyzzzz, g_y_0_xxxy_xyzzzzz, g_y_0_xxxy_xzzzzz, g_y_0_xxxy_xzzzzzz, g_y_0_xxxy_yyyyyy, g_y_0_xxxy_yyyyyz, g_y_0_xxxy_yyyyzz, g_y_0_xxxy_yyyzzz, g_y_0_xxxy_yyzzzz, g_y_0_xxxy_yzzzzz, g_y_0_xxxy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxy_xxxxxx[k] = -g_y_0_xxxy_xxxxxx[k] * cd_x[k] + g_y_0_xxxy_xxxxxxx[k];

                g_y_0_xxxxy_xxxxxy[k] = -g_y_0_xxxy_xxxxxy[k] * cd_x[k] + g_y_0_xxxy_xxxxxxy[k];

                g_y_0_xxxxy_xxxxxz[k] = -g_y_0_xxxy_xxxxxz[k] * cd_x[k] + g_y_0_xxxy_xxxxxxz[k];

                g_y_0_xxxxy_xxxxyy[k] = -g_y_0_xxxy_xxxxyy[k] * cd_x[k] + g_y_0_xxxy_xxxxxyy[k];

                g_y_0_xxxxy_xxxxyz[k] = -g_y_0_xxxy_xxxxyz[k] * cd_x[k] + g_y_0_xxxy_xxxxxyz[k];

                g_y_0_xxxxy_xxxxzz[k] = -g_y_0_xxxy_xxxxzz[k] * cd_x[k] + g_y_0_xxxy_xxxxxzz[k];

                g_y_0_xxxxy_xxxyyy[k] = -g_y_0_xxxy_xxxyyy[k] * cd_x[k] + g_y_0_xxxy_xxxxyyy[k];

                g_y_0_xxxxy_xxxyyz[k] = -g_y_0_xxxy_xxxyyz[k] * cd_x[k] + g_y_0_xxxy_xxxxyyz[k];

                g_y_0_xxxxy_xxxyzz[k] = -g_y_0_xxxy_xxxyzz[k] * cd_x[k] + g_y_0_xxxy_xxxxyzz[k];

                g_y_0_xxxxy_xxxzzz[k] = -g_y_0_xxxy_xxxzzz[k] * cd_x[k] + g_y_0_xxxy_xxxxzzz[k];

                g_y_0_xxxxy_xxyyyy[k] = -g_y_0_xxxy_xxyyyy[k] * cd_x[k] + g_y_0_xxxy_xxxyyyy[k];

                g_y_0_xxxxy_xxyyyz[k] = -g_y_0_xxxy_xxyyyz[k] * cd_x[k] + g_y_0_xxxy_xxxyyyz[k];

                g_y_0_xxxxy_xxyyzz[k] = -g_y_0_xxxy_xxyyzz[k] * cd_x[k] + g_y_0_xxxy_xxxyyzz[k];

                g_y_0_xxxxy_xxyzzz[k] = -g_y_0_xxxy_xxyzzz[k] * cd_x[k] + g_y_0_xxxy_xxxyzzz[k];

                g_y_0_xxxxy_xxzzzz[k] = -g_y_0_xxxy_xxzzzz[k] * cd_x[k] + g_y_0_xxxy_xxxzzzz[k];

                g_y_0_xxxxy_xyyyyy[k] = -g_y_0_xxxy_xyyyyy[k] * cd_x[k] + g_y_0_xxxy_xxyyyyy[k];

                g_y_0_xxxxy_xyyyyz[k] = -g_y_0_xxxy_xyyyyz[k] * cd_x[k] + g_y_0_xxxy_xxyyyyz[k];

                g_y_0_xxxxy_xyyyzz[k] = -g_y_0_xxxy_xyyyzz[k] * cd_x[k] + g_y_0_xxxy_xxyyyzz[k];

                g_y_0_xxxxy_xyyzzz[k] = -g_y_0_xxxy_xyyzzz[k] * cd_x[k] + g_y_0_xxxy_xxyyzzz[k];

                g_y_0_xxxxy_xyzzzz[k] = -g_y_0_xxxy_xyzzzz[k] * cd_x[k] + g_y_0_xxxy_xxyzzzz[k];

                g_y_0_xxxxy_xzzzzz[k] = -g_y_0_xxxy_xzzzzz[k] * cd_x[k] + g_y_0_xxxy_xxzzzzz[k];

                g_y_0_xxxxy_yyyyyy[k] = -g_y_0_xxxy_yyyyyy[k] * cd_x[k] + g_y_0_xxxy_xyyyyyy[k];

                g_y_0_xxxxy_yyyyyz[k] = -g_y_0_xxxy_yyyyyz[k] * cd_x[k] + g_y_0_xxxy_xyyyyyz[k];

                g_y_0_xxxxy_yyyyzz[k] = -g_y_0_xxxy_yyyyzz[k] * cd_x[k] + g_y_0_xxxy_xyyyyzz[k];

                g_y_0_xxxxy_yyyzzz[k] = -g_y_0_xxxy_yyyzzz[k] * cd_x[k] + g_y_0_xxxy_xyyyzzz[k];

                g_y_0_xxxxy_yyzzzz[k] = -g_y_0_xxxy_yyzzzz[k] * cd_x[k] + g_y_0_xxxy_xyyzzzz[k];

                g_y_0_xxxxy_yzzzzz[k] = -g_y_0_xxxy_yzzzzz[k] * cd_x[k] + g_y_0_xxxy_xyzzzzz[k];

                g_y_0_xxxxy_zzzzzz[k] = -g_y_0_xxxy_zzzzzz[k] * cd_x[k] + g_y_0_xxxy_xzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 56);

            auto g_y_0_xxxxz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 57);

            auto g_y_0_xxxxz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 58);

            auto g_y_0_xxxxz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 59);

            auto g_y_0_xxxxz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 60);

            auto g_y_0_xxxxz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 61);

            auto g_y_0_xxxxz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 62);

            auto g_y_0_xxxxz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 63);

            auto g_y_0_xxxxz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 64);

            auto g_y_0_xxxxz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 65);

            auto g_y_0_xxxxz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 66);

            auto g_y_0_xxxxz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 67);

            auto g_y_0_xxxxz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 68);

            auto g_y_0_xxxxz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 69);

            auto g_y_0_xxxxz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 70);

            auto g_y_0_xxxxz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 71);

            auto g_y_0_xxxxz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 72);

            auto g_y_0_xxxxz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 73);

            auto g_y_0_xxxxz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 74);

            auto g_y_0_xxxxz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 75);

            auto g_y_0_xxxxz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 76);

            auto g_y_0_xxxxz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 77);

            auto g_y_0_xxxxz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 78);

            auto g_y_0_xxxxz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 79);

            auto g_y_0_xxxxz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 80);

            auto g_y_0_xxxxz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 81);

            auto g_y_0_xxxxz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 82);

            auto g_y_0_xxxxz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxz_xxxxxx, g_y_0_xxxxz_xxxxxy, g_y_0_xxxxz_xxxxxz, g_y_0_xxxxz_xxxxyy, g_y_0_xxxxz_xxxxyz, g_y_0_xxxxz_xxxxzz, g_y_0_xxxxz_xxxyyy, g_y_0_xxxxz_xxxyyz, g_y_0_xxxxz_xxxyzz, g_y_0_xxxxz_xxxzzz, g_y_0_xxxxz_xxyyyy, g_y_0_xxxxz_xxyyyz, g_y_0_xxxxz_xxyyzz, g_y_0_xxxxz_xxyzzz, g_y_0_xxxxz_xxzzzz, g_y_0_xxxxz_xyyyyy, g_y_0_xxxxz_xyyyyz, g_y_0_xxxxz_xyyyzz, g_y_0_xxxxz_xyyzzz, g_y_0_xxxxz_xyzzzz, g_y_0_xxxxz_xzzzzz, g_y_0_xxxxz_yyyyyy, g_y_0_xxxxz_yyyyyz, g_y_0_xxxxz_yyyyzz, g_y_0_xxxxz_yyyzzz, g_y_0_xxxxz_yyzzzz, g_y_0_xxxxz_yzzzzz, g_y_0_xxxxz_zzzzzz, g_y_0_xxxz_xxxxxx, g_y_0_xxxz_xxxxxxx, g_y_0_xxxz_xxxxxxy, g_y_0_xxxz_xxxxxxz, g_y_0_xxxz_xxxxxy, g_y_0_xxxz_xxxxxyy, g_y_0_xxxz_xxxxxyz, g_y_0_xxxz_xxxxxz, g_y_0_xxxz_xxxxxzz, g_y_0_xxxz_xxxxyy, g_y_0_xxxz_xxxxyyy, g_y_0_xxxz_xxxxyyz, g_y_0_xxxz_xxxxyz, g_y_0_xxxz_xxxxyzz, g_y_0_xxxz_xxxxzz, g_y_0_xxxz_xxxxzzz, g_y_0_xxxz_xxxyyy, g_y_0_xxxz_xxxyyyy, g_y_0_xxxz_xxxyyyz, g_y_0_xxxz_xxxyyz, g_y_0_xxxz_xxxyyzz, g_y_0_xxxz_xxxyzz, g_y_0_xxxz_xxxyzzz, g_y_0_xxxz_xxxzzz, g_y_0_xxxz_xxxzzzz, g_y_0_xxxz_xxyyyy, g_y_0_xxxz_xxyyyyy, g_y_0_xxxz_xxyyyyz, g_y_0_xxxz_xxyyyz, g_y_0_xxxz_xxyyyzz, g_y_0_xxxz_xxyyzz, g_y_0_xxxz_xxyyzzz, g_y_0_xxxz_xxyzzz, g_y_0_xxxz_xxyzzzz, g_y_0_xxxz_xxzzzz, g_y_0_xxxz_xxzzzzz, g_y_0_xxxz_xyyyyy, g_y_0_xxxz_xyyyyyy, g_y_0_xxxz_xyyyyyz, g_y_0_xxxz_xyyyyz, g_y_0_xxxz_xyyyyzz, g_y_0_xxxz_xyyyzz, g_y_0_xxxz_xyyyzzz, g_y_0_xxxz_xyyzzz, g_y_0_xxxz_xyyzzzz, g_y_0_xxxz_xyzzzz, g_y_0_xxxz_xyzzzzz, g_y_0_xxxz_xzzzzz, g_y_0_xxxz_xzzzzzz, g_y_0_xxxz_yyyyyy, g_y_0_xxxz_yyyyyz, g_y_0_xxxz_yyyyzz, g_y_0_xxxz_yyyzzz, g_y_0_xxxz_yyzzzz, g_y_0_xxxz_yzzzzz, g_y_0_xxxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxz_xxxxxx[k] = -g_y_0_xxxz_xxxxxx[k] * cd_x[k] + g_y_0_xxxz_xxxxxxx[k];

                g_y_0_xxxxz_xxxxxy[k] = -g_y_0_xxxz_xxxxxy[k] * cd_x[k] + g_y_0_xxxz_xxxxxxy[k];

                g_y_0_xxxxz_xxxxxz[k] = -g_y_0_xxxz_xxxxxz[k] * cd_x[k] + g_y_0_xxxz_xxxxxxz[k];

                g_y_0_xxxxz_xxxxyy[k] = -g_y_0_xxxz_xxxxyy[k] * cd_x[k] + g_y_0_xxxz_xxxxxyy[k];

                g_y_0_xxxxz_xxxxyz[k] = -g_y_0_xxxz_xxxxyz[k] * cd_x[k] + g_y_0_xxxz_xxxxxyz[k];

                g_y_0_xxxxz_xxxxzz[k] = -g_y_0_xxxz_xxxxzz[k] * cd_x[k] + g_y_0_xxxz_xxxxxzz[k];

                g_y_0_xxxxz_xxxyyy[k] = -g_y_0_xxxz_xxxyyy[k] * cd_x[k] + g_y_0_xxxz_xxxxyyy[k];

                g_y_0_xxxxz_xxxyyz[k] = -g_y_0_xxxz_xxxyyz[k] * cd_x[k] + g_y_0_xxxz_xxxxyyz[k];

                g_y_0_xxxxz_xxxyzz[k] = -g_y_0_xxxz_xxxyzz[k] * cd_x[k] + g_y_0_xxxz_xxxxyzz[k];

                g_y_0_xxxxz_xxxzzz[k] = -g_y_0_xxxz_xxxzzz[k] * cd_x[k] + g_y_0_xxxz_xxxxzzz[k];

                g_y_0_xxxxz_xxyyyy[k] = -g_y_0_xxxz_xxyyyy[k] * cd_x[k] + g_y_0_xxxz_xxxyyyy[k];

                g_y_0_xxxxz_xxyyyz[k] = -g_y_0_xxxz_xxyyyz[k] * cd_x[k] + g_y_0_xxxz_xxxyyyz[k];

                g_y_0_xxxxz_xxyyzz[k] = -g_y_0_xxxz_xxyyzz[k] * cd_x[k] + g_y_0_xxxz_xxxyyzz[k];

                g_y_0_xxxxz_xxyzzz[k] = -g_y_0_xxxz_xxyzzz[k] * cd_x[k] + g_y_0_xxxz_xxxyzzz[k];

                g_y_0_xxxxz_xxzzzz[k] = -g_y_0_xxxz_xxzzzz[k] * cd_x[k] + g_y_0_xxxz_xxxzzzz[k];

                g_y_0_xxxxz_xyyyyy[k] = -g_y_0_xxxz_xyyyyy[k] * cd_x[k] + g_y_0_xxxz_xxyyyyy[k];

                g_y_0_xxxxz_xyyyyz[k] = -g_y_0_xxxz_xyyyyz[k] * cd_x[k] + g_y_0_xxxz_xxyyyyz[k];

                g_y_0_xxxxz_xyyyzz[k] = -g_y_0_xxxz_xyyyzz[k] * cd_x[k] + g_y_0_xxxz_xxyyyzz[k];

                g_y_0_xxxxz_xyyzzz[k] = -g_y_0_xxxz_xyyzzz[k] * cd_x[k] + g_y_0_xxxz_xxyyzzz[k];

                g_y_0_xxxxz_xyzzzz[k] = -g_y_0_xxxz_xyzzzz[k] * cd_x[k] + g_y_0_xxxz_xxyzzzz[k];

                g_y_0_xxxxz_xzzzzz[k] = -g_y_0_xxxz_xzzzzz[k] * cd_x[k] + g_y_0_xxxz_xxzzzzz[k];

                g_y_0_xxxxz_yyyyyy[k] = -g_y_0_xxxz_yyyyyy[k] * cd_x[k] + g_y_0_xxxz_xyyyyyy[k];

                g_y_0_xxxxz_yyyyyz[k] = -g_y_0_xxxz_yyyyyz[k] * cd_x[k] + g_y_0_xxxz_xyyyyyz[k];

                g_y_0_xxxxz_yyyyzz[k] = -g_y_0_xxxz_yyyyzz[k] * cd_x[k] + g_y_0_xxxz_xyyyyzz[k];

                g_y_0_xxxxz_yyyzzz[k] = -g_y_0_xxxz_yyyzzz[k] * cd_x[k] + g_y_0_xxxz_xyyyzzz[k];

                g_y_0_xxxxz_yyzzzz[k] = -g_y_0_xxxz_yyzzzz[k] * cd_x[k] + g_y_0_xxxz_xyyzzzz[k];

                g_y_0_xxxxz_yzzzzz[k] = -g_y_0_xxxz_yzzzzz[k] * cd_x[k] + g_y_0_xxxz_xyzzzzz[k];

                g_y_0_xxxxz_zzzzzz[k] = -g_y_0_xxxz_zzzzzz[k] * cd_x[k] + g_y_0_xxxz_xzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyy_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 84);

            auto g_y_0_xxxyy_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 85);

            auto g_y_0_xxxyy_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 86);

            auto g_y_0_xxxyy_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 87);

            auto g_y_0_xxxyy_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 88);

            auto g_y_0_xxxyy_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 89);

            auto g_y_0_xxxyy_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 90);

            auto g_y_0_xxxyy_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 91);

            auto g_y_0_xxxyy_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 92);

            auto g_y_0_xxxyy_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 93);

            auto g_y_0_xxxyy_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 94);

            auto g_y_0_xxxyy_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 95);

            auto g_y_0_xxxyy_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 96);

            auto g_y_0_xxxyy_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 97);

            auto g_y_0_xxxyy_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 98);

            auto g_y_0_xxxyy_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 99);

            auto g_y_0_xxxyy_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 100);

            auto g_y_0_xxxyy_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 101);

            auto g_y_0_xxxyy_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 102);

            auto g_y_0_xxxyy_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 103);

            auto g_y_0_xxxyy_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 104);

            auto g_y_0_xxxyy_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 105);

            auto g_y_0_xxxyy_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 106);

            auto g_y_0_xxxyy_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 107);

            auto g_y_0_xxxyy_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 108);

            auto g_y_0_xxxyy_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 109);

            auto g_y_0_xxxyy_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 110);

            auto g_y_0_xxxyy_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 111);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyy_xxxxxx, g_y_0_xxxyy_xxxxxy, g_y_0_xxxyy_xxxxxz, g_y_0_xxxyy_xxxxyy, g_y_0_xxxyy_xxxxyz, g_y_0_xxxyy_xxxxzz, g_y_0_xxxyy_xxxyyy, g_y_0_xxxyy_xxxyyz, g_y_0_xxxyy_xxxyzz, g_y_0_xxxyy_xxxzzz, g_y_0_xxxyy_xxyyyy, g_y_0_xxxyy_xxyyyz, g_y_0_xxxyy_xxyyzz, g_y_0_xxxyy_xxyzzz, g_y_0_xxxyy_xxzzzz, g_y_0_xxxyy_xyyyyy, g_y_0_xxxyy_xyyyyz, g_y_0_xxxyy_xyyyzz, g_y_0_xxxyy_xyyzzz, g_y_0_xxxyy_xyzzzz, g_y_0_xxxyy_xzzzzz, g_y_0_xxxyy_yyyyyy, g_y_0_xxxyy_yyyyyz, g_y_0_xxxyy_yyyyzz, g_y_0_xxxyy_yyyzzz, g_y_0_xxxyy_yyzzzz, g_y_0_xxxyy_yzzzzz, g_y_0_xxxyy_zzzzzz, g_y_0_xxyy_xxxxxx, g_y_0_xxyy_xxxxxxx, g_y_0_xxyy_xxxxxxy, g_y_0_xxyy_xxxxxxz, g_y_0_xxyy_xxxxxy, g_y_0_xxyy_xxxxxyy, g_y_0_xxyy_xxxxxyz, g_y_0_xxyy_xxxxxz, g_y_0_xxyy_xxxxxzz, g_y_0_xxyy_xxxxyy, g_y_0_xxyy_xxxxyyy, g_y_0_xxyy_xxxxyyz, g_y_0_xxyy_xxxxyz, g_y_0_xxyy_xxxxyzz, g_y_0_xxyy_xxxxzz, g_y_0_xxyy_xxxxzzz, g_y_0_xxyy_xxxyyy, g_y_0_xxyy_xxxyyyy, g_y_0_xxyy_xxxyyyz, g_y_0_xxyy_xxxyyz, g_y_0_xxyy_xxxyyzz, g_y_0_xxyy_xxxyzz, g_y_0_xxyy_xxxyzzz, g_y_0_xxyy_xxxzzz, g_y_0_xxyy_xxxzzzz, g_y_0_xxyy_xxyyyy, g_y_0_xxyy_xxyyyyy, g_y_0_xxyy_xxyyyyz, g_y_0_xxyy_xxyyyz, g_y_0_xxyy_xxyyyzz, g_y_0_xxyy_xxyyzz, g_y_0_xxyy_xxyyzzz, g_y_0_xxyy_xxyzzz, g_y_0_xxyy_xxyzzzz, g_y_0_xxyy_xxzzzz, g_y_0_xxyy_xxzzzzz, g_y_0_xxyy_xyyyyy, g_y_0_xxyy_xyyyyyy, g_y_0_xxyy_xyyyyyz, g_y_0_xxyy_xyyyyz, g_y_0_xxyy_xyyyyzz, g_y_0_xxyy_xyyyzz, g_y_0_xxyy_xyyyzzz, g_y_0_xxyy_xyyzzz, g_y_0_xxyy_xyyzzzz, g_y_0_xxyy_xyzzzz, g_y_0_xxyy_xyzzzzz, g_y_0_xxyy_xzzzzz, g_y_0_xxyy_xzzzzzz, g_y_0_xxyy_yyyyyy, g_y_0_xxyy_yyyyyz, g_y_0_xxyy_yyyyzz, g_y_0_xxyy_yyyzzz, g_y_0_xxyy_yyzzzz, g_y_0_xxyy_yzzzzz, g_y_0_xxyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyy_xxxxxx[k] = -g_y_0_xxyy_xxxxxx[k] * cd_x[k] + g_y_0_xxyy_xxxxxxx[k];

                g_y_0_xxxyy_xxxxxy[k] = -g_y_0_xxyy_xxxxxy[k] * cd_x[k] + g_y_0_xxyy_xxxxxxy[k];

                g_y_0_xxxyy_xxxxxz[k] = -g_y_0_xxyy_xxxxxz[k] * cd_x[k] + g_y_0_xxyy_xxxxxxz[k];

                g_y_0_xxxyy_xxxxyy[k] = -g_y_0_xxyy_xxxxyy[k] * cd_x[k] + g_y_0_xxyy_xxxxxyy[k];

                g_y_0_xxxyy_xxxxyz[k] = -g_y_0_xxyy_xxxxyz[k] * cd_x[k] + g_y_0_xxyy_xxxxxyz[k];

                g_y_0_xxxyy_xxxxzz[k] = -g_y_0_xxyy_xxxxzz[k] * cd_x[k] + g_y_0_xxyy_xxxxxzz[k];

                g_y_0_xxxyy_xxxyyy[k] = -g_y_0_xxyy_xxxyyy[k] * cd_x[k] + g_y_0_xxyy_xxxxyyy[k];

                g_y_0_xxxyy_xxxyyz[k] = -g_y_0_xxyy_xxxyyz[k] * cd_x[k] + g_y_0_xxyy_xxxxyyz[k];

                g_y_0_xxxyy_xxxyzz[k] = -g_y_0_xxyy_xxxyzz[k] * cd_x[k] + g_y_0_xxyy_xxxxyzz[k];

                g_y_0_xxxyy_xxxzzz[k] = -g_y_0_xxyy_xxxzzz[k] * cd_x[k] + g_y_0_xxyy_xxxxzzz[k];

                g_y_0_xxxyy_xxyyyy[k] = -g_y_0_xxyy_xxyyyy[k] * cd_x[k] + g_y_0_xxyy_xxxyyyy[k];

                g_y_0_xxxyy_xxyyyz[k] = -g_y_0_xxyy_xxyyyz[k] * cd_x[k] + g_y_0_xxyy_xxxyyyz[k];

                g_y_0_xxxyy_xxyyzz[k] = -g_y_0_xxyy_xxyyzz[k] * cd_x[k] + g_y_0_xxyy_xxxyyzz[k];

                g_y_0_xxxyy_xxyzzz[k] = -g_y_0_xxyy_xxyzzz[k] * cd_x[k] + g_y_0_xxyy_xxxyzzz[k];

                g_y_0_xxxyy_xxzzzz[k] = -g_y_0_xxyy_xxzzzz[k] * cd_x[k] + g_y_0_xxyy_xxxzzzz[k];

                g_y_0_xxxyy_xyyyyy[k] = -g_y_0_xxyy_xyyyyy[k] * cd_x[k] + g_y_0_xxyy_xxyyyyy[k];

                g_y_0_xxxyy_xyyyyz[k] = -g_y_0_xxyy_xyyyyz[k] * cd_x[k] + g_y_0_xxyy_xxyyyyz[k];

                g_y_0_xxxyy_xyyyzz[k] = -g_y_0_xxyy_xyyyzz[k] * cd_x[k] + g_y_0_xxyy_xxyyyzz[k];

                g_y_0_xxxyy_xyyzzz[k] = -g_y_0_xxyy_xyyzzz[k] * cd_x[k] + g_y_0_xxyy_xxyyzzz[k];

                g_y_0_xxxyy_xyzzzz[k] = -g_y_0_xxyy_xyzzzz[k] * cd_x[k] + g_y_0_xxyy_xxyzzzz[k];

                g_y_0_xxxyy_xzzzzz[k] = -g_y_0_xxyy_xzzzzz[k] * cd_x[k] + g_y_0_xxyy_xxzzzzz[k];

                g_y_0_xxxyy_yyyyyy[k] = -g_y_0_xxyy_yyyyyy[k] * cd_x[k] + g_y_0_xxyy_xyyyyyy[k];

                g_y_0_xxxyy_yyyyyz[k] = -g_y_0_xxyy_yyyyyz[k] * cd_x[k] + g_y_0_xxyy_xyyyyyz[k];

                g_y_0_xxxyy_yyyyzz[k] = -g_y_0_xxyy_yyyyzz[k] * cd_x[k] + g_y_0_xxyy_xyyyyzz[k];

                g_y_0_xxxyy_yyyzzz[k] = -g_y_0_xxyy_yyyzzz[k] * cd_x[k] + g_y_0_xxyy_xyyyzzz[k];

                g_y_0_xxxyy_yyzzzz[k] = -g_y_0_xxyy_yyzzzz[k] * cd_x[k] + g_y_0_xxyy_xyyzzzz[k];

                g_y_0_xxxyy_yzzzzz[k] = -g_y_0_xxyy_yzzzzz[k] * cd_x[k] + g_y_0_xxyy_xyzzzzz[k];

                g_y_0_xxxyy_zzzzzz[k] = -g_y_0_xxyy_zzzzzz[k] * cd_x[k] + g_y_0_xxyy_xzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 112);

            auto g_y_0_xxxyz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 113);

            auto g_y_0_xxxyz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 114);

            auto g_y_0_xxxyz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 115);

            auto g_y_0_xxxyz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 116);

            auto g_y_0_xxxyz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 117);

            auto g_y_0_xxxyz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 118);

            auto g_y_0_xxxyz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 119);

            auto g_y_0_xxxyz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 120);

            auto g_y_0_xxxyz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 121);

            auto g_y_0_xxxyz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 122);

            auto g_y_0_xxxyz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 123);

            auto g_y_0_xxxyz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 124);

            auto g_y_0_xxxyz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 125);

            auto g_y_0_xxxyz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 126);

            auto g_y_0_xxxyz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 127);

            auto g_y_0_xxxyz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 128);

            auto g_y_0_xxxyz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 129);

            auto g_y_0_xxxyz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 130);

            auto g_y_0_xxxyz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 131);

            auto g_y_0_xxxyz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 132);

            auto g_y_0_xxxyz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 133);

            auto g_y_0_xxxyz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 134);

            auto g_y_0_xxxyz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 135);

            auto g_y_0_xxxyz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 136);

            auto g_y_0_xxxyz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 137);

            auto g_y_0_xxxyz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 138);

            auto g_y_0_xxxyz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 139);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyz_xxxxxx, g_y_0_xxxyz_xxxxxy, g_y_0_xxxyz_xxxxxz, g_y_0_xxxyz_xxxxyy, g_y_0_xxxyz_xxxxyz, g_y_0_xxxyz_xxxxzz, g_y_0_xxxyz_xxxyyy, g_y_0_xxxyz_xxxyyz, g_y_0_xxxyz_xxxyzz, g_y_0_xxxyz_xxxzzz, g_y_0_xxxyz_xxyyyy, g_y_0_xxxyz_xxyyyz, g_y_0_xxxyz_xxyyzz, g_y_0_xxxyz_xxyzzz, g_y_0_xxxyz_xxzzzz, g_y_0_xxxyz_xyyyyy, g_y_0_xxxyz_xyyyyz, g_y_0_xxxyz_xyyyzz, g_y_0_xxxyz_xyyzzz, g_y_0_xxxyz_xyzzzz, g_y_0_xxxyz_xzzzzz, g_y_0_xxxyz_yyyyyy, g_y_0_xxxyz_yyyyyz, g_y_0_xxxyz_yyyyzz, g_y_0_xxxyz_yyyzzz, g_y_0_xxxyz_yyzzzz, g_y_0_xxxyz_yzzzzz, g_y_0_xxxyz_zzzzzz, g_y_0_xxyz_xxxxxx, g_y_0_xxyz_xxxxxxx, g_y_0_xxyz_xxxxxxy, g_y_0_xxyz_xxxxxxz, g_y_0_xxyz_xxxxxy, g_y_0_xxyz_xxxxxyy, g_y_0_xxyz_xxxxxyz, g_y_0_xxyz_xxxxxz, g_y_0_xxyz_xxxxxzz, g_y_0_xxyz_xxxxyy, g_y_0_xxyz_xxxxyyy, g_y_0_xxyz_xxxxyyz, g_y_0_xxyz_xxxxyz, g_y_0_xxyz_xxxxyzz, g_y_0_xxyz_xxxxzz, g_y_0_xxyz_xxxxzzz, g_y_0_xxyz_xxxyyy, g_y_0_xxyz_xxxyyyy, g_y_0_xxyz_xxxyyyz, g_y_0_xxyz_xxxyyz, g_y_0_xxyz_xxxyyzz, g_y_0_xxyz_xxxyzz, g_y_0_xxyz_xxxyzzz, g_y_0_xxyz_xxxzzz, g_y_0_xxyz_xxxzzzz, g_y_0_xxyz_xxyyyy, g_y_0_xxyz_xxyyyyy, g_y_0_xxyz_xxyyyyz, g_y_0_xxyz_xxyyyz, g_y_0_xxyz_xxyyyzz, g_y_0_xxyz_xxyyzz, g_y_0_xxyz_xxyyzzz, g_y_0_xxyz_xxyzzz, g_y_0_xxyz_xxyzzzz, g_y_0_xxyz_xxzzzz, g_y_0_xxyz_xxzzzzz, g_y_0_xxyz_xyyyyy, g_y_0_xxyz_xyyyyyy, g_y_0_xxyz_xyyyyyz, g_y_0_xxyz_xyyyyz, g_y_0_xxyz_xyyyyzz, g_y_0_xxyz_xyyyzz, g_y_0_xxyz_xyyyzzz, g_y_0_xxyz_xyyzzz, g_y_0_xxyz_xyyzzzz, g_y_0_xxyz_xyzzzz, g_y_0_xxyz_xyzzzzz, g_y_0_xxyz_xzzzzz, g_y_0_xxyz_xzzzzzz, g_y_0_xxyz_yyyyyy, g_y_0_xxyz_yyyyyz, g_y_0_xxyz_yyyyzz, g_y_0_xxyz_yyyzzz, g_y_0_xxyz_yyzzzz, g_y_0_xxyz_yzzzzz, g_y_0_xxyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyz_xxxxxx[k] = -g_y_0_xxyz_xxxxxx[k] * cd_x[k] + g_y_0_xxyz_xxxxxxx[k];

                g_y_0_xxxyz_xxxxxy[k] = -g_y_0_xxyz_xxxxxy[k] * cd_x[k] + g_y_0_xxyz_xxxxxxy[k];

                g_y_0_xxxyz_xxxxxz[k] = -g_y_0_xxyz_xxxxxz[k] * cd_x[k] + g_y_0_xxyz_xxxxxxz[k];

                g_y_0_xxxyz_xxxxyy[k] = -g_y_0_xxyz_xxxxyy[k] * cd_x[k] + g_y_0_xxyz_xxxxxyy[k];

                g_y_0_xxxyz_xxxxyz[k] = -g_y_0_xxyz_xxxxyz[k] * cd_x[k] + g_y_0_xxyz_xxxxxyz[k];

                g_y_0_xxxyz_xxxxzz[k] = -g_y_0_xxyz_xxxxzz[k] * cd_x[k] + g_y_0_xxyz_xxxxxzz[k];

                g_y_0_xxxyz_xxxyyy[k] = -g_y_0_xxyz_xxxyyy[k] * cd_x[k] + g_y_0_xxyz_xxxxyyy[k];

                g_y_0_xxxyz_xxxyyz[k] = -g_y_0_xxyz_xxxyyz[k] * cd_x[k] + g_y_0_xxyz_xxxxyyz[k];

                g_y_0_xxxyz_xxxyzz[k] = -g_y_0_xxyz_xxxyzz[k] * cd_x[k] + g_y_0_xxyz_xxxxyzz[k];

                g_y_0_xxxyz_xxxzzz[k] = -g_y_0_xxyz_xxxzzz[k] * cd_x[k] + g_y_0_xxyz_xxxxzzz[k];

                g_y_0_xxxyz_xxyyyy[k] = -g_y_0_xxyz_xxyyyy[k] * cd_x[k] + g_y_0_xxyz_xxxyyyy[k];

                g_y_0_xxxyz_xxyyyz[k] = -g_y_0_xxyz_xxyyyz[k] * cd_x[k] + g_y_0_xxyz_xxxyyyz[k];

                g_y_0_xxxyz_xxyyzz[k] = -g_y_0_xxyz_xxyyzz[k] * cd_x[k] + g_y_0_xxyz_xxxyyzz[k];

                g_y_0_xxxyz_xxyzzz[k] = -g_y_0_xxyz_xxyzzz[k] * cd_x[k] + g_y_0_xxyz_xxxyzzz[k];

                g_y_0_xxxyz_xxzzzz[k] = -g_y_0_xxyz_xxzzzz[k] * cd_x[k] + g_y_0_xxyz_xxxzzzz[k];

                g_y_0_xxxyz_xyyyyy[k] = -g_y_0_xxyz_xyyyyy[k] * cd_x[k] + g_y_0_xxyz_xxyyyyy[k];

                g_y_0_xxxyz_xyyyyz[k] = -g_y_0_xxyz_xyyyyz[k] * cd_x[k] + g_y_0_xxyz_xxyyyyz[k];

                g_y_0_xxxyz_xyyyzz[k] = -g_y_0_xxyz_xyyyzz[k] * cd_x[k] + g_y_0_xxyz_xxyyyzz[k];

                g_y_0_xxxyz_xyyzzz[k] = -g_y_0_xxyz_xyyzzz[k] * cd_x[k] + g_y_0_xxyz_xxyyzzz[k];

                g_y_0_xxxyz_xyzzzz[k] = -g_y_0_xxyz_xyzzzz[k] * cd_x[k] + g_y_0_xxyz_xxyzzzz[k];

                g_y_0_xxxyz_xzzzzz[k] = -g_y_0_xxyz_xzzzzz[k] * cd_x[k] + g_y_0_xxyz_xxzzzzz[k];

                g_y_0_xxxyz_yyyyyy[k] = -g_y_0_xxyz_yyyyyy[k] * cd_x[k] + g_y_0_xxyz_xyyyyyy[k];

                g_y_0_xxxyz_yyyyyz[k] = -g_y_0_xxyz_yyyyyz[k] * cd_x[k] + g_y_0_xxyz_xyyyyyz[k];

                g_y_0_xxxyz_yyyyzz[k] = -g_y_0_xxyz_yyyyzz[k] * cd_x[k] + g_y_0_xxyz_xyyyyzz[k];

                g_y_0_xxxyz_yyyzzz[k] = -g_y_0_xxyz_yyyzzz[k] * cd_x[k] + g_y_0_xxyz_xyyyzzz[k];

                g_y_0_xxxyz_yyzzzz[k] = -g_y_0_xxyz_yyzzzz[k] * cd_x[k] + g_y_0_xxyz_xyyzzzz[k];

                g_y_0_xxxyz_yzzzzz[k] = -g_y_0_xxyz_yzzzzz[k] * cd_x[k] + g_y_0_xxyz_xyzzzzz[k];

                g_y_0_xxxyz_zzzzzz[k] = -g_y_0_xxyz_zzzzzz[k] * cd_x[k] + g_y_0_xxyz_xzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 140);

            auto g_y_0_xxxzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 141);

            auto g_y_0_xxxzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 142);

            auto g_y_0_xxxzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 143);

            auto g_y_0_xxxzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 144);

            auto g_y_0_xxxzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 145);

            auto g_y_0_xxxzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 146);

            auto g_y_0_xxxzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 147);

            auto g_y_0_xxxzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 148);

            auto g_y_0_xxxzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 149);

            auto g_y_0_xxxzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 150);

            auto g_y_0_xxxzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 151);

            auto g_y_0_xxxzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 152);

            auto g_y_0_xxxzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 153);

            auto g_y_0_xxxzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 154);

            auto g_y_0_xxxzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 155);

            auto g_y_0_xxxzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 156);

            auto g_y_0_xxxzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 157);

            auto g_y_0_xxxzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 158);

            auto g_y_0_xxxzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 159);

            auto g_y_0_xxxzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 160);

            auto g_y_0_xxxzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 161);

            auto g_y_0_xxxzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 162);

            auto g_y_0_xxxzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 163);

            auto g_y_0_xxxzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 164);

            auto g_y_0_xxxzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 165);

            auto g_y_0_xxxzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 166);

            auto g_y_0_xxxzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_x, g_y_0_xxxzz_xxxxxx, g_y_0_xxxzz_xxxxxy, g_y_0_xxxzz_xxxxxz, g_y_0_xxxzz_xxxxyy, g_y_0_xxxzz_xxxxyz, g_y_0_xxxzz_xxxxzz, g_y_0_xxxzz_xxxyyy, g_y_0_xxxzz_xxxyyz, g_y_0_xxxzz_xxxyzz, g_y_0_xxxzz_xxxzzz, g_y_0_xxxzz_xxyyyy, g_y_0_xxxzz_xxyyyz, g_y_0_xxxzz_xxyyzz, g_y_0_xxxzz_xxyzzz, g_y_0_xxxzz_xxzzzz, g_y_0_xxxzz_xyyyyy, g_y_0_xxxzz_xyyyyz, g_y_0_xxxzz_xyyyzz, g_y_0_xxxzz_xyyzzz, g_y_0_xxxzz_xyzzzz, g_y_0_xxxzz_xzzzzz, g_y_0_xxxzz_yyyyyy, g_y_0_xxxzz_yyyyyz, g_y_0_xxxzz_yyyyzz, g_y_0_xxxzz_yyyzzz, g_y_0_xxxzz_yyzzzz, g_y_0_xxxzz_yzzzzz, g_y_0_xxxzz_zzzzzz, g_y_0_xxzz_xxxxxx, g_y_0_xxzz_xxxxxxx, g_y_0_xxzz_xxxxxxy, g_y_0_xxzz_xxxxxxz, g_y_0_xxzz_xxxxxy, g_y_0_xxzz_xxxxxyy, g_y_0_xxzz_xxxxxyz, g_y_0_xxzz_xxxxxz, g_y_0_xxzz_xxxxxzz, g_y_0_xxzz_xxxxyy, g_y_0_xxzz_xxxxyyy, g_y_0_xxzz_xxxxyyz, g_y_0_xxzz_xxxxyz, g_y_0_xxzz_xxxxyzz, g_y_0_xxzz_xxxxzz, g_y_0_xxzz_xxxxzzz, g_y_0_xxzz_xxxyyy, g_y_0_xxzz_xxxyyyy, g_y_0_xxzz_xxxyyyz, g_y_0_xxzz_xxxyyz, g_y_0_xxzz_xxxyyzz, g_y_0_xxzz_xxxyzz, g_y_0_xxzz_xxxyzzz, g_y_0_xxzz_xxxzzz, g_y_0_xxzz_xxxzzzz, g_y_0_xxzz_xxyyyy, g_y_0_xxzz_xxyyyyy, g_y_0_xxzz_xxyyyyz, g_y_0_xxzz_xxyyyz, g_y_0_xxzz_xxyyyzz, g_y_0_xxzz_xxyyzz, g_y_0_xxzz_xxyyzzz, g_y_0_xxzz_xxyzzz, g_y_0_xxzz_xxyzzzz, g_y_0_xxzz_xxzzzz, g_y_0_xxzz_xxzzzzz, g_y_0_xxzz_xyyyyy, g_y_0_xxzz_xyyyyyy, g_y_0_xxzz_xyyyyyz, g_y_0_xxzz_xyyyyz, g_y_0_xxzz_xyyyyzz, g_y_0_xxzz_xyyyzz, g_y_0_xxzz_xyyyzzz, g_y_0_xxzz_xyyzzz, g_y_0_xxzz_xyyzzzz, g_y_0_xxzz_xyzzzz, g_y_0_xxzz_xyzzzzz, g_y_0_xxzz_xzzzzz, g_y_0_xxzz_xzzzzzz, g_y_0_xxzz_yyyyyy, g_y_0_xxzz_yyyyyz, g_y_0_xxzz_yyyyzz, g_y_0_xxzz_yyyzzz, g_y_0_xxzz_yyzzzz, g_y_0_xxzz_yzzzzz, g_y_0_xxzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzz_xxxxxx[k] = -g_y_0_xxzz_xxxxxx[k] * cd_x[k] + g_y_0_xxzz_xxxxxxx[k];

                g_y_0_xxxzz_xxxxxy[k] = -g_y_0_xxzz_xxxxxy[k] * cd_x[k] + g_y_0_xxzz_xxxxxxy[k];

                g_y_0_xxxzz_xxxxxz[k] = -g_y_0_xxzz_xxxxxz[k] * cd_x[k] + g_y_0_xxzz_xxxxxxz[k];

                g_y_0_xxxzz_xxxxyy[k] = -g_y_0_xxzz_xxxxyy[k] * cd_x[k] + g_y_0_xxzz_xxxxxyy[k];

                g_y_0_xxxzz_xxxxyz[k] = -g_y_0_xxzz_xxxxyz[k] * cd_x[k] + g_y_0_xxzz_xxxxxyz[k];

                g_y_0_xxxzz_xxxxzz[k] = -g_y_0_xxzz_xxxxzz[k] * cd_x[k] + g_y_0_xxzz_xxxxxzz[k];

                g_y_0_xxxzz_xxxyyy[k] = -g_y_0_xxzz_xxxyyy[k] * cd_x[k] + g_y_0_xxzz_xxxxyyy[k];

                g_y_0_xxxzz_xxxyyz[k] = -g_y_0_xxzz_xxxyyz[k] * cd_x[k] + g_y_0_xxzz_xxxxyyz[k];

                g_y_0_xxxzz_xxxyzz[k] = -g_y_0_xxzz_xxxyzz[k] * cd_x[k] + g_y_0_xxzz_xxxxyzz[k];

                g_y_0_xxxzz_xxxzzz[k] = -g_y_0_xxzz_xxxzzz[k] * cd_x[k] + g_y_0_xxzz_xxxxzzz[k];

                g_y_0_xxxzz_xxyyyy[k] = -g_y_0_xxzz_xxyyyy[k] * cd_x[k] + g_y_0_xxzz_xxxyyyy[k];

                g_y_0_xxxzz_xxyyyz[k] = -g_y_0_xxzz_xxyyyz[k] * cd_x[k] + g_y_0_xxzz_xxxyyyz[k];

                g_y_0_xxxzz_xxyyzz[k] = -g_y_0_xxzz_xxyyzz[k] * cd_x[k] + g_y_0_xxzz_xxxyyzz[k];

                g_y_0_xxxzz_xxyzzz[k] = -g_y_0_xxzz_xxyzzz[k] * cd_x[k] + g_y_0_xxzz_xxxyzzz[k];

                g_y_0_xxxzz_xxzzzz[k] = -g_y_0_xxzz_xxzzzz[k] * cd_x[k] + g_y_0_xxzz_xxxzzzz[k];

                g_y_0_xxxzz_xyyyyy[k] = -g_y_0_xxzz_xyyyyy[k] * cd_x[k] + g_y_0_xxzz_xxyyyyy[k];

                g_y_0_xxxzz_xyyyyz[k] = -g_y_0_xxzz_xyyyyz[k] * cd_x[k] + g_y_0_xxzz_xxyyyyz[k];

                g_y_0_xxxzz_xyyyzz[k] = -g_y_0_xxzz_xyyyzz[k] * cd_x[k] + g_y_0_xxzz_xxyyyzz[k];

                g_y_0_xxxzz_xyyzzz[k] = -g_y_0_xxzz_xyyzzz[k] * cd_x[k] + g_y_0_xxzz_xxyyzzz[k];

                g_y_0_xxxzz_xyzzzz[k] = -g_y_0_xxzz_xyzzzz[k] * cd_x[k] + g_y_0_xxzz_xxyzzzz[k];

                g_y_0_xxxzz_xzzzzz[k] = -g_y_0_xxzz_xzzzzz[k] * cd_x[k] + g_y_0_xxzz_xxzzzzz[k];

                g_y_0_xxxzz_yyyyyy[k] = -g_y_0_xxzz_yyyyyy[k] * cd_x[k] + g_y_0_xxzz_xyyyyyy[k];

                g_y_0_xxxzz_yyyyyz[k] = -g_y_0_xxzz_yyyyyz[k] * cd_x[k] + g_y_0_xxzz_xyyyyyz[k];

                g_y_0_xxxzz_yyyyzz[k] = -g_y_0_xxzz_yyyyzz[k] * cd_x[k] + g_y_0_xxzz_xyyyyzz[k];

                g_y_0_xxxzz_yyyzzz[k] = -g_y_0_xxzz_yyyzzz[k] * cd_x[k] + g_y_0_xxzz_xyyyzzz[k];

                g_y_0_xxxzz_yyzzzz[k] = -g_y_0_xxzz_yyzzzz[k] * cd_x[k] + g_y_0_xxzz_xyyzzzz[k];

                g_y_0_xxxzz_yzzzzz[k] = -g_y_0_xxzz_yzzzzz[k] * cd_x[k] + g_y_0_xxzz_xyzzzzz[k];

                g_y_0_xxxzz_zzzzzz[k] = -g_y_0_xxzz_zzzzzz[k] * cd_x[k] + g_y_0_xxzz_xzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 168);

            auto g_y_0_xxyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 169);

            auto g_y_0_xxyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 170);

            auto g_y_0_xxyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 171);

            auto g_y_0_xxyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 172);

            auto g_y_0_xxyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 173);

            auto g_y_0_xxyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 174);

            auto g_y_0_xxyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 175);

            auto g_y_0_xxyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 176);

            auto g_y_0_xxyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 177);

            auto g_y_0_xxyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 178);

            auto g_y_0_xxyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 179);

            auto g_y_0_xxyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 180);

            auto g_y_0_xxyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 181);

            auto g_y_0_xxyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 182);

            auto g_y_0_xxyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 183);

            auto g_y_0_xxyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 184);

            auto g_y_0_xxyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 185);

            auto g_y_0_xxyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 186);

            auto g_y_0_xxyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 187);

            auto g_y_0_xxyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 188);

            auto g_y_0_xxyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 189);

            auto g_y_0_xxyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 190);

            auto g_y_0_xxyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 191);

            auto g_y_0_xxyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 192);

            auto g_y_0_xxyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 193);

            auto g_y_0_xxyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 194);

            auto g_y_0_xxyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 195);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyy_xxxxxx, g_y_0_xxyyy_xxxxxy, g_y_0_xxyyy_xxxxxz, g_y_0_xxyyy_xxxxyy, g_y_0_xxyyy_xxxxyz, g_y_0_xxyyy_xxxxzz, g_y_0_xxyyy_xxxyyy, g_y_0_xxyyy_xxxyyz, g_y_0_xxyyy_xxxyzz, g_y_0_xxyyy_xxxzzz, g_y_0_xxyyy_xxyyyy, g_y_0_xxyyy_xxyyyz, g_y_0_xxyyy_xxyyzz, g_y_0_xxyyy_xxyzzz, g_y_0_xxyyy_xxzzzz, g_y_0_xxyyy_xyyyyy, g_y_0_xxyyy_xyyyyz, g_y_0_xxyyy_xyyyzz, g_y_0_xxyyy_xyyzzz, g_y_0_xxyyy_xyzzzz, g_y_0_xxyyy_xzzzzz, g_y_0_xxyyy_yyyyyy, g_y_0_xxyyy_yyyyyz, g_y_0_xxyyy_yyyyzz, g_y_0_xxyyy_yyyzzz, g_y_0_xxyyy_yyzzzz, g_y_0_xxyyy_yzzzzz, g_y_0_xxyyy_zzzzzz, g_y_0_xyyy_xxxxxx, g_y_0_xyyy_xxxxxxx, g_y_0_xyyy_xxxxxxy, g_y_0_xyyy_xxxxxxz, g_y_0_xyyy_xxxxxy, g_y_0_xyyy_xxxxxyy, g_y_0_xyyy_xxxxxyz, g_y_0_xyyy_xxxxxz, g_y_0_xyyy_xxxxxzz, g_y_0_xyyy_xxxxyy, g_y_0_xyyy_xxxxyyy, g_y_0_xyyy_xxxxyyz, g_y_0_xyyy_xxxxyz, g_y_0_xyyy_xxxxyzz, g_y_0_xyyy_xxxxzz, g_y_0_xyyy_xxxxzzz, g_y_0_xyyy_xxxyyy, g_y_0_xyyy_xxxyyyy, g_y_0_xyyy_xxxyyyz, g_y_0_xyyy_xxxyyz, g_y_0_xyyy_xxxyyzz, g_y_0_xyyy_xxxyzz, g_y_0_xyyy_xxxyzzz, g_y_0_xyyy_xxxzzz, g_y_0_xyyy_xxxzzzz, g_y_0_xyyy_xxyyyy, g_y_0_xyyy_xxyyyyy, g_y_0_xyyy_xxyyyyz, g_y_0_xyyy_xxyyyz, g_y_0_xyyy_xxyyyzz, g_y_0_xyyy_xxyyzz, g_y_0_xyyy_xxyyzzz, g_y_0_xyyy_xxyzzz, g_y_0_xyyy_xxyzzzz, g_y_0_xyyy_xxzzzz, g_y_0_xyyy_xxzzzzz, g_y_0_xyyy_xyyyyy, g_y_0_xyyy_xyyyyyy, g_y_0_xyyy_xyyyyyz, g_y_0_xyyy_xyyyyz, g_y_0_xyyy_xyyyyzz, g_y_0_xyyy_xyyyzz, g_y_0_xyyy_xyyyzzz, g_y_0_xyyy_xyyzzz, g_y_0_xyyy_xyyzzzz, g_y_0_xyyy_xyzzzz, g_y_0_xyyy_xyzzzzz, g_y_0_xyyy_xzzzzz, g_y_0_xyyy_xzzzzzz, g_y_0_xyyy_yyyyyy, g_y_0_xyyy_yyyyyz, g_y_0_xyyy_yyyyzz, g_y_0_xyyy_yyyzzz, g_y_0_xyyy_yyzzzz, g_y_0_xyyy_yzzzzz, g_y_0_xyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyy_xxxxxx[k] = -g_y_0_xyyy_xxxxxx[k] * cd_x[k] + g_y_0_xyyy_xxxxxxx[k];

                g_y_0_xxyyy_xxxxxy[k] = -g_y_0_xyyy_xxxxxy[k] * cd_x[k] + g_y_0_xyyy_xxxxxxy[k];

                g_y_0_xxyyy_xxxxxz[k] = -g_y_0_xyyy_xxxxxz[k] * cd_x[k] + g_y_0_xyyy_xxxxxxz[k];

                g_y_0_xxyyy_xxxxyy[k] = -g_y_0_xyyy_xxxxyy[k] * cd_x[k] + g_y_0_xyyy_xxxxxyy[k];

                g_y_0_xxyyy_xxxxyz[k] = -g_y_0_xyyy_xxxxyz[k] * cd_x[k] + g_y_0_xyyy_xxxxxyz[k];

                g_y_0_xxyyy_xxxxzz[k] = -g_y_0_xyyy_xxxxzz[k] * cd_x[k] + g_y_0_xyyy_xxxxxzz[k];

                g_y_0_xxyyy_xxxyyy[k] = -g_y_0_xyyy_xxxyyy[k] * cd_x[k] + g_y_0_xyyy_xxxxyyy[k];

                g_y_0_xxyyy_xxxyyz[k] = -g_y_0_xyyy_xxxyyz[k] * cd_x[k] + g_y_0_xyyy_xxxxyyz[k];

                g_y_0_xxyyy_xxxyzz[k] = -g_y_0_xyyy_xxxyzz[k] * cd_x[k] + g_y_0_xyyy_xxxxyzz[k];

                g_y_0_xxyyy_xxxzzz[k] = -g_y_0_xyyy_xxxzzz[k] * cd_x[k] + g_y_0_xyyy_xxxxzzz[k];

                g_y_0_xxyyy_xxyyyy[k] = -g_y_0_xyyy_xxyyyy[k] * cd_x[k] + g_y_0_xyyy_xxxyyyy[k];

                g_y_0_xxyyy_xxyyyz[k] = -g_y_0_xyyy_xxyyyz[k] * cd_x[k] + g_y_0_xyyy_xxxyyyz[k];

                g_y_0_xxyyy_xxyyzz[k] = -g_y_0_xyyy_xxyyzz[k] * cd_x[k] + g_y_0_xyyy_xxxyyzz[k];

                g_y_0_xxyyy_xxyzzz[k] = -g_y_0_xyyy_xxyzzz[k] * cd_x[k] + g_y_0_xyyy_xxxyzzz[k];

                g_y_0_xxyyy_xxzzzz[k] = -g_y_0_xyyy_xxzzzz[k] * cd_x[k] + g_y_0_xyyy_xxxzzzz[k];

                g_y_0_xxyyy_xyyyyy[k] = -g_y_0_xyyy_xyyyyy[k] * cd_x[k] + g_y_0_xyyy_xxyyyyy[k];

                g_y_0_xxyyy_xyyyyz[k] = -g_y_0_xyyy_xyyyyz[k] * cd_x[k] + g_y_0_xyyy_xxyyyyz[k];

                g_y_0_xxyyy_xyyyzz[k] = -g_y_0_xyyy_xyyyzz[k] * cd_x[k] + g_y_0_xyyy_xxyyyzz[k];

                g_y_0_xxyyy_xyyzzz[k] = -g_y_0_xyyy_xyyzzz[k] * cd_x[k] + g_y_0_xyyy_xxyyzzz[k];

                g_y_0_xxyyy_xyzzzz[k] = -g_y_0_xyyy_xyzzzz[k] * cd_x[k] + g_y_0_xyyy_xxyzzzz[k];

                g_y_0_xxyyy_xzzzzz[k] = -g_y_0_xyyy_xzzzzz[k] * cd_x[k] + g_y_0_xyyy_xxzzzzz[k];

                g_y_0_xxyyy_yyyyyy[k] = -g_y_0_xyyy_yyyyyy[k] * cd_x[k] + g_y_0_xyyy_xyyyyyy[k];

                g_y_0_xxyyy_yyyyyz[k] = -g_y_0_xyyy_yyyyyz[k] * cd_x[k] + g_y_0_xyyy_xyyyyyz[k];

                g_y_0_xxyyy_yyyyzz[k] = -g_y_0_xyyy_yyyyzz[k] * cd_x[k] + g_y_0_xyyy_xyyyyzz[k];

                g_y_0_xxyyy_yyyzzz[k] = -g_y_0_xyyy_yyyzzz[k] * cd_x[k] + g_y_0_xyyy_xyyyzzz[k];

                g_y_0_xxyyy_yyzzzz[k] = -g_y_0_xyyy_yyzzzz[k] * cd_x[k] + g_y_0_xyyy_xyyzzzz[k];

                g_y_0_xxyyy_yzzzzz[k] = -g_y_0_xyyy_yzzzzz[k] * cd_x[k] + g_y_0_xyyy_xyzzzzz[k];

                g_y_0_xxyyy_zzzzzz[k] = -g_y_0_xyyy_zzzzzz[k] * cd_x[k] + g_y_0_xyyy_xzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 196);

            auto g_y_0_xxyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 197);

            auto g_y_0_xxyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 198);

            auto g_y_0_xxyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 199);

            auto g_y_0_xxyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 200);

            auto g_y_0_xxyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 201);

            auto g_y_0_xxyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 202);

            auto g_y_0_xxyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 203);

            auto g_y_0_xxyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 204);

            auto g_y_0_xxyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 205);

            auto g_y_0_xxyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 206);

            auto g_y_0_xxyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 207);

            auto g_y_0_xxyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 208);

            auto g_y_0_xxyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 209);

            auto g_y_0_xxyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 210);

            auto g_y_0_xxyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 211);

            auto g_y_0_xxyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 212);

            auto g_y_0_xxyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 213);

            auto g_y_0_xxyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 214);

            auto g_y_0_xxyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 215);

            auto g_y_0_xxyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 216);

            auto g_y_0_xxyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 217);

            auto g_y_0_xxyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 218);

            auto g_y_0_xxyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 219);

            auto g_y_0_xxyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 220);

            auto g_y_0_xxyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 221);

            auto g_y_0_xxyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 222);

            auto g_y_0_xxyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 223);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyz_xxxxxx, g_y_0_xxyyz_xxxxxy, g_y_0_xxyyz_xxxxxz, g_y_0_xxyyz_xxxxyy, g_y_0_xxyyz_xxxxyz, g_y_0_xxyyz_xxxxzz, g_y_0_xxyyz_xxxyyy, g_y_0_xxyyz_xxxyyz, g_y_0_xxyyz_xxxyzz, g_y_0_xxyyz_xxxzzz, g_y_0_xxyyz_xxyyyy, g_y_0_xxyyz_xxyyyz, g_y_0_xxyyz_xxyyzz, g_y_0_xxyyz_xxyzzz, g_y_0_xxyyz_xxzzzz, g_y_0_xxyyz_xyyyyy, g_y_0_xxyyz_xyyyyz, g_y_0_xxyyz_xyyyzz, g_y_0_xxyyz_xyyzzz, g_y_0_xxyyz_xyzzzz, g_y_0_xxyyz_xzzzzz, g_y_0_xxyyz_yyyyyy, g_y_0_xxyyz_yyyyyz, g_y_0_xxyyz_yyyyzz, g_y_0_xxyyz_yyyzzz, g_y_0_xxyyz_yyzzzz, g_y_0_xxyyz_yzzzzz, g_y_0_xxyyz_zzzzzz, g_y_0_xyyz_xxxxxx, g_y_0_xyyz_xxxxxxx, g_y_0_xyyz_xxxxxxy, g_y_0_xyyz_xxxxxxz, g_y_0_xyyz_xxxxxy, g_y_0_xyyz_xxxxxyy, g_y_0_xyyz_xxxxxyz, g_y_0_xyyz_xxxxxz, g_y_0_xyyz_xxxxxzz, g_y_0_xyyz_xxxxyy, g_y_0_xyyz_xxxxyyy, g_y_0_xyyz_xxxxyyz, g_y_0_xyyz_xxxxyz, g_y_0_xyyz_xxxxyzz, g_y_0_xyyz_xxxxzz, g_y_0_xyyz_xxxxzzz, g_y_0_xyyz_xxxyyy, g_y_0_xyyz_xxxyyyy, g_y_0_xyyz_xxxyyyz, g_y_0_xyyz_xxxyyz, g_y_0_xyyz_xxxyyzz, g_y_0_xyyz_xxxyzz, g_y_0_xyyz_xxxyzzz, g_y_0_xyyz_xxxzzz, g_y_0_xyyz_xxxzzzz, g_y_0_xyyz_xxyyyy, g_y_0_xyyz_xxyyyyy, g_y_0_xyyz_xxyyyyz, g_y_0_xyyz_xxyyyz, g_y_0_xyyz_xxyyyzz, g_y_0_xyyz_xxyyzz, g_y_0_xyyz_xxyyzzz, g_y_0_xyyz_xxyzzz, g_y_0_xyyz_xxyzzzz, g_y_0_xyyz_xxzzzz, g_y_0_xyyz_xxzzzzz, g_y_0_xyyz_xyyyyy, g_y_0_xyyz_xyyyyyy, g_y_0_xyyz_xyyyyyz, g_y_0_xyyz_xyyyyz, g_y_0_xyyz_xyyyyzz, g_y_0_xyyz_xyyyzz, g_y_0_xyyz_xyyyzzz, g_y_0_xyyz_xyyzzz, g_y_0_xyyz_xyyzzzz, g_y_0_xyyz_xyzzzz, g_y_0_xyyz_xyzzzzz, g_y_0_xyyz_xzzzzz, g_y_0_xyyz_xzzzzzz, g_y_0_xyyz_yyyyyy, g_y_0_xyyz_yyyyyz, g_y_0_xyyz_yyyyzz, g_y_0_xyyz_yyyzzz, g_y_0_xyyz_yyzzzz, g_y_0_xyyz_yzzzzz, g_y_0_xyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyz_xxxxxx[k] = -g_y_0_xyyz_xxxxxx[k] * cd_x[k] + g_y_0_xyyz_xxxxxxx[k];

                g_y_0_xxyyz_xxxxxy[k] = -g_y_0_xyyz_xxxxxy[k] * cd_x[k] + g_y_0_xyyz_xxxxxxy[k];

                g_y_0_xxyyz_xxxxxz[k] = -g_y_0_xyyz_xxxxxz[k] * cd_x[k] + g_y_0_xyyz_xxxxxxz[k];

                g_y_0_xxyyz_xxxxyy[k] = -g_y_0_xyyz_xxxxyy[k] * cd_x[k] + g_y_0_xyyz_xxxxxyy[k];

                g_y_0_xxyyz_xxxxyz[k] = -g_y_0_xyyz_xxxxyz[k] * cd_x[k] + g_y_0_xyyz_xxxxxyz[k];

                g_y_0_xxyyz_xxxxzz[k] = -g_y_0_xyyz_xxxxzz[k] * cd_x[k] + g_y_0_xyyz_xxxxxzz[k];

                g_y_0_xxyyz_xxxyyy[k] = -g_y_0_xyyz_xxxyyy[k] * cd_x[k] + g_y_0_xyyz_xxxxyyy[k];

                g_y_0_xxyyz_xxxyyz[k] = -g_y_0_xyyz_xxxyyz[k] * cd_x[k] + g_y_0_xyyz_xxxxyyz[k];

                g_y_0_xxyyz_xxxyzz[k] = -g_y_0_xyyz_xxxyzz[k] * cd_x[k] + g_y_0_xyyz_xxxxyzz[k];

                g_y_0_xxyyz_xxxzzz[k] = -g_y_0_xyyz_xxxzzz[k] * cd_x[k] + g_y_0_xyyz_xxxxzzz[k];

                g_y_0_xxyyz_xxyyyy[k] = -g_y_0_xyyz_xxyyyy[k] * cd_x[k] + g_y_0_xyyz_xxxyyyy[k];

                g_y_0_xxyyz_xxyyyz[k] = -g_y_0_xyyz_xxyyyz[k] * cd_x[k] + g_y_0_xyyz_xxxyyyz[k];

                g_y_0_xxyyz_xxyyzz[k] = -g_y_0_xyyz_xxyyzz[k] * cd_x[k] + g_y_0_xyyz_xxxyyzz[k];

                g_y_0_xxyyz_xxyzzz[k] = -g_y_0_xyyz_xxyzzz[k] * cd_x[k] + g_y_0_xyyz_xxxyzzz[k];

                g_y_0_xxyyz_xxzzzz[k] = -g_y_0_xyyz_xxzzzz[k] * cd_x[k] + g_y_0_xyyz_xxxzzzz[k];

                g_y_0_xxyyz_xyyyyy[k] = -g_y_0_xyyz_xyyyyy[k] * cd_x[k] + g_y_0_xyyz_xxyyyyy[k];

                g_y_0_xxyyz_xyyyyz[k] = -g_y_0_xyyz_xyyyyz[k] * cd_x[k] + g_y_0_xyyz_xxyyyyz[k];

                g_y_0_xxyyz_xyyyzz[k] = -g_y_0_xyyz_xyyyzz[k] * cd_x[k] + g_y_0_xyyz_xxyyyzz[k];

                g_y_0_xxyyz_xyyzzz[k] = -g_y_0_xyyz_xyyzzz[k] * cd_x[k] + g_y_0_xyyz_xxyyzzz[k];

                g_y_0_xxyyz_xyzzzz[k] = -g_y_0_xyyz_xyzzzz[k] * cd_x[k] + g_y_0_xyyz_xxyzzzz[k];

                g_y_0_xxyyz_xzzzzz[k] = -g_y_0_xyyz_xzzzzz[k] * cd_x[k] + g_y_0_xyyz_xxzzzzz[k];

                g_y_0_xxyyz_yyyyyy[k] = -g_y_0_xyyz_yyyyyy[k] * cd_x[k] + g_y_0_xyyz_xyyyyyy[k];

                g_y_0_xxyyz_yyyyyz[k] = -g_y_0_xyyz_yyyyyz[k] * cd_x[k] + g_y_0_xyyz_xyyyyyz[k];

                g_y_0_xxyyz_yyyyzz[k] = -g_y_0_xyyz_yyyyzz[k] * cd_x[k] + g_y_0_xyyz_xyyyyzz[k];

                g_y_0_xxyyz_yyyzzz[k] = -g_y_0_xyyz_yyyzzz[k] * cd_x[k] + g_y_0_xyyz_xyyyzzz[k];

                g_y_0_xxyyz_yyzzzz[k] = -g_y_0_xyyz_yyzzzz[k] * cd_x[k] + g_y_0_xyyz_xyyzzzz[k];

                g_y_0_xxyyz_yzzzzz[k] = -g_y_0_xyyz_yzzzzz[k] * cd_x[k] + g_y_0_xyyz_xyzzzzz[k];

                g_y_0_xxyyz_zzzzzz[k] = -g_y_0_xyyz_zzzzzz[k] * cd_x[k] + g_y_0_xyyz_xzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 224);

            auto g_y_0_xxyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 225);

            auto g_y_0_xxyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 226);

            auto g_y_0_xxyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 227);

            auto g_y_0_xxyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 228);

            auto g_y_0_xxyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 229);

            auto g_y_0_xxyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 230);

            auto g_y_0_xxyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 231);

            auto g_y_0_xxyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 232);

            auto g_y_0_xxyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 233);

            auto g_y_0_xxyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 234);

            auto g_y_0_xxyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 235);

            auto g_y_0_xxyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 236);

            auto g_y_0_xxyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 237);

            auto g_y_0_xxyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 238);

            auto g_y_0_xxyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 239);

            auto g_y_0_xxyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 240);

            auto g_y_0_xxyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 241);

            auto g_y_0_xxyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 242);

            auto g_y_0_xxyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 243);

            auto g_y_0_xxyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 244);

            auto g_y_0_xxyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 245);

            auto g_y_0_xxyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 246);

            auto g_y_0_xxyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 247);

            auto g_y_0_xxyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 248);

            auto g_y_0_xxyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 249);

            auto g_y_0_xxyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 250);

            auto g_y_0_xxyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 251);

            #pragma omp simd aligned(cd_x, g_y_0_xxyzz_xxxxxx, g_y_0_xxyzz_xxxxxy, g_y_0_xxyzz_xxxxxz, g_y_0_xxyzz_xxxxyy, g_y_0_xxyzz_xxxxyz, g_y_0_xxyzz_xxxxzz, g_y_0_xxyzz_xxxyyy, g_y_0_xxyzz_xxxyyz, g_y_0_xxyzz_xxxyzz, g_y_0_xxyzz_xxxzzz, g_y_0_xxyzz_xxyyyy, g_y_0_xxyzz_xxyyyz, g_y_0_xxyzz_xxyyzz, g_y_0_xxyzz_xxyzzz, g_y_0_xxyzz_xxzzzz, g_y_0_xxyzz_xyyyyy, g_y_0_xxyzz_xyyyyz, g_y_0_xxyzz_xyyyzz, g_y_0_xxyzz_xyyzzz, g_y_0_xxyzz_xyzzzz, g_y_0_xxyzz_xzzzzz, g_y_0_xxyzz_yyyyyy, g_y_0_xxyzz_yyyyyz, g_y_0_xxyzz_yyyyzz, g_y_0_xxyzz_yyyzzz, g_y_0_xxyzz_yyzzzz, g_y_0_xxyzz_yzzzzz, g_y_0_xxyzz_zzzzzz, g_y_0_xyzz_xxxxxx, g_y_0_xyzz_xxxxxxx, g_y_0_xyzz_xxxxxxy, g_y_0_xyzz_xxxxxxz, g_y_0_xyzz_xxxxxy, g_y_0_xyzz_xxxxxyy, g_y_0_xyzz_xxxxxyz, g_y_0_xyzz_xxxxxz, g_y_0_xyzz_xxxxxzz, g_y_0_xyzz_xxxxyy, g_y_0_xyzz_xxxxyyy, g_y_0_xyzz_xxxxyyz, g_y_0_xyzz_xxxxyz, g_y_0_xyzz_xxxxyzz, g_y_0_xyzz_xxxxzz, g_y_0_xyzz_xxxxzzz, g_y_0_xyzz_xxxyyy, g_y_0_xyzz_xxxyyyy, g_y_0_xyzz_xxxyyyz, g_y_0_xyzz_xxxyyz, g_y_0_xyzz_xxxyyzz, g_y_0_xyzz_xxxyzz, g_y_0_xyzz_xxxyzzz, g_y_0_xyzz_xxxzzz, g_y_0_xyzz_xxxzzzz, g_y_0_xyzz_xxyyyy, g_y_0_xyzz_xxyyyyy, g_y_0_xyzz_xxyyyyz, g_y_0_xyzz_xxyyyz, g_y_0_xyzz_xxyyyzz, g_y_0_xyzz_xxyyzz, g_y_0_xyzz_xxyyzzz, g_y_0_xyzz_xxyzzz, g_y_0_xyzz_xxyzzzz, g_y_0_xyzz_xxzzzz, g_y_0_xyzz_xxzzzzz, g_y_0_xyzz_xyyyyy, g_y_0_xyzz_xyyyyyy, g_y_0_xyzz_xyyyyyz, g_y_0_xyzz_xyyyyz, g_y_0_xyzz_xyyyyzz, g_y_0_xyzz_xyyyzz, g_y_0_xyzz_xyyyzzz, g_y_0_xyzz_xyyzzz, g_y_0_xyzz_xyyzzzz, g_y_0_xyzz_xyzzzz, g_y_0_xyzz_xyzzzzz, g_y_0_xyzz_xzzzzz, g_y_0_xyzz_xzzzzzz, g_y_0_xyzz_yyyyyy, g_y_0_xyzz_yyyyyz, g_y_0_xyzz_yyyyzz, g_y_0_xyzz_yyyzzz, g_y_0_xyzz_yyzzzz, g_y_0_xyzz_yzzzzz, g_y_0_xyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzz_xxxxxx[k] = -g_y_0_xyzz_xxxxxx[k] * cd_x[k] + g_y_0_xyzz_xxxxxxx[k];

                g_y_0_xxyzz_xxxxxy[k] = -g_y_0_xyzz_xxxxxy[k] * cd_x[k] + g_y_0_xyzz_xxxxxxy[k];

                g_y_0_xxyzz_xxxxxz[k] = -g_y_0_xyzz_xxxxxz[k] * cd_x[k] + g_y_0_xyzz_xxxxxxz[k];

                g_y_0_xxyzz_xxxxyy[k] = -g_y_0_xyzz_xxxxyy[k] * cd_x[k] + g_y_0_xyzz_xxxxxyy[k];

                g_y_0_xxyzz_xxxxyz[k] = -g_y_0_xyzz_xxxxyz[k] * cd_x[k] + g_y_0_xyzz_xxxxxyz[k];

                g_y_0_xxyzz_xxxxzz[k] = -g_y_0_xyzz_xxxxzz[k] * cd_x[k] + g_y_0_xyzz_xxxxxzz[k];

                g_y_0_xxyzz_xxxyyy[k] = -g_y_0_xyzz_xxxyyy[k] * cd_x[k] + g_y_0_xyzz_xxxxyyy[k];

                g_y_0_xxyzz_xxxyyz[k] = -g_y_0_xyzz_xxxyyz[k] * cd_x[k] + g_y_0_xyzz_xxxxyyz[k];

                g_y_0_xxyzz_xxxyzz[k] = -g_y_0_xyzz_xxxyzz[k] * cd_x[k] + g_y_0_xyzz_xxxxyzz[k];

                g_y_0_xxyzz_xxxzzz[k] = -g_y_0_xyzz_xxxzzz[k] * cd_x[k] + g_y_0_xyzz_xxxxzzz[k];

                g_y_0_xxyzz_xxyyyy[k] = -g_y_0_xyzz_xxyyyy[k] * cd_x[k] + g_y_0_xyzz_xxxyyyy[k];

                g_y_0_xxyzz_xxyyyz[k] = -g_y_0_xyzz_xxyyyz[k] * cd_x[k] + g_y_0_xyzz_xxxyyyz[k];

                g_y_0_xxyzz_xxyyzz[k] = -g_y_0_xyzz_xxyyzz[k] * cd_x[k] + g_y_0_xyzz_xxxyyzz[k];

                g_y_0_xxyzz_xxyzzz[k] = -g_y_0_xyzz_xxyzzz[k] * cd_x[k] + g_y_0_xyzz_xxxyzzz[k];

                g_y_0_xxyzz_xxzzzz[k] = -g_y_0_xyzz_xxzzzz[k] * cd_x[k] + g_y_0_xyzz_xxxzzzz[k];

                g_y_0_xxyzz_xyyyyy[k] = -g_y_0_xyzz_xyyyyy[k] * cd_x[k] + g_y_0_xyzz_xxyyyyy[k];

                g_y_0_xxyzz_xyyyyz[k] = -g_y_0_xyzz_xyyyyz[k] * cd_x[k] + g_y_0_xyzz_xxyyyyz[k];

                g_y_0_xxyzz_xyyyzz[k] = -g_y_0_xyzz_xyyyzz[k] * cd_x[k] + g_y_0_xyzz_xxyyyzz[k];

                g_y_0_xxyzz_xyyzzz[k] = -g_y_0_xyzz_xyyzzz[k] * cd_x[k] + g_y_0_xyzz_xxyyzzz[k];

                g_y_0_xxyzz_xyzzzz[k] = -g_y_0_xyzz_xyzzzz[k] * cd_x[k] + g_y_0_xyzz_xxyzzzz[k];

                g_y_0_xxyzz_xzzzzz[k] = -g_y_0_xyzz_xzzzzz[k] * cd_x[k] + g_y_0_xyzz_xxzzzzz[k];

                g_y_0_xxyzz_yyyyyy[k] = -g_y_0_xyzz_yyyyyy[k] * cd_x[k] + g_y_0_xyzz_xyyyyyy[k];

                g_y_0_xxyzz_yyyyyz[k] = -g_y_0_xyzz_yyyyyz[k] * cd_x[k] + g_y_0_xyzz_xyyyyyz[k];

                g_y_0_xxyzz_yyyyzz[k] = -g_y_0_xyzz_yyyyzz[k] * cd_x[k] + g_y_0_xyzz_xyyyyzz[k];

                g_y_0_xxyzz_yyyzzz[k] = -g_y_0_xyzz_yyyzzz[k] * cd_x[k] + g_y_0_xyzz_xyyyzzz[k];

                g_y_0_xxyzz_yyzzzz[k] = -g_y_0_xyzz_yyzzzz[k] * cd_x[k] + g_y_0_xyzz_xyyzzzz[k];

                g_y_0_xxyzz_yzzzzz[k] = -g_y_0_xyzz_yzzzzz[k] * cd_x[k] + g_y_0_xyzz_xyzzzzz[k];

                g_y_0_xxyzz_zzzzzz[k] = -g_y_0_xyzz_zzzzzz[k] * cd_x[k] + g_y_0_xyzz_xzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 252);

            auto g_y_0_xxzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 253);

            auto g_y_0_xxzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 254);

            auto g_y_0_xxzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 255);

            auto g_y_0_xxzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 256);

            auto g_y_0_xxzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 257);

            auto g_y_0_xxzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 258);

            auto g_y_0_xxzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 259);

            auto g_y_0_xxzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 260);

            auto g_y_0_xxzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 261);

            auto g_y_0_xxzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 262);

            auto g_y_0_xxzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 263);

            auto g_y_0_xxzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 264);

            auto g_y_0_xxzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 265);

            auto g_y_0_xxzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 266);

            auto g_y_0_xxzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 267);

            auto g_y_0_xxzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 268);

            auto g_y_0_xxzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 269);

            auto g_y_0_xxzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 270);

            auto g_y_0_xxzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 271);

            auto g_y_0_xxzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 272);

            auto g_y_0_xxzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 273);

            auto g_y_0_xxzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 274);

            auto g_y_0_xxzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 275);

            auto g_y_0_xxzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 276);

            auto g_y_0_xxzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 277);

            auto g_y_0_xxzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 278);

            auto g_y_0_xxzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 279);

            #pragma omp simd aligned(cd_x, g_y_0_xxzzz_xxxxxx, g_y_0_xxzzz_xxxxxy, g_y_0_xxzzz_xxxxxz, g_y_0_xxzzz_xxxxyy, g_y_0_xxzzz_xxxxyz, g_y_0_xxzzz_xxxxzz, g_y_0_xxzzz_xxxyyy, g_y_0_xxzzz_xxxyyz, g_y_0_xxzzz_xxxyzz, g_y_0_xxzzz_xxxzzz, g_y_0_xxzzz_xxyyyy, g_y_0_xxzzz_xxyyyz, g_y_0_xxzzz_xxyyzz, g_y_0_xxzzz_xxyzzz, g_y_0_xxzzz_xxzzzz, g_y_0_xxzzz_xyyyyy, g_y_0_xxzzz_xyyyyz, g_y_0_xxzzz_xyyyzz, g_y_0_xxzzz_xyyzzz, g_y_0_xxzzz_xyzzzz, g_y_0_xxzzz_xzzzzz, g_y_0_xxzzz_yyyyyy, g_y_0_xxzzz_yyyyyz, g_y_0_xxzzz_yyyyzz, g_y_0_xxzzz_yyyzzz, g_y_0_xxzzz_yyzzzz, g_y_0_xxzzz_yzzzzz, g_y_0_xxzzz_zzzzzz, g_y_0_xzzz_xxxxxx, g_y_0_xzzz_xxxxxxx, g_y_0_xzzz_xxxxxxy, g_y_0_xzzz_xxxxxxz, g_y_0_xzzz_xxxxxy, g_y_0_xzzz_xxxxxyy, g_y_0_xzzz_xxxxxyz, g_y_0_xzzz_xxxxxz, g_y_0_xzzz_xxxxxzz, g_y_0_xzzz_xxxxyy, g_y_0_xzzz_xxxxyyy, g_y_0_xzzz_xxxxyyz, g_y_0_xzzz_xxxxyz, g_y_0_xzzz_xxxxyzz, g_y_0_xzzz_xxxxzz, g_y_0_xzzz_xxxxzzz, g_y_0_xzzz_xxxyyy, g_y_0_xzzz_xxxyyyy, g_y_0_xzzz_xxxyyyz, g_y_0_xzzz_xxxyyz, g_y_0_xzzz_xxxyyzz, g_y_0_xzzz_xxxyzz, g_y_0_xzzz_xxxyzzz, g_y_0_xzzz_xxxzzz, g_y_0_xzzz_xxxzzzz, g_y_0_xzzz_xxyyyy, g_y_0_xzzz_xxyyyyy, g_y_0_xzzz_xxyyyyz, g_y_0_xzzz_xxyyyz, g_y_0_xzzz_xxyyyzz, g_y_0_xzzz_xxyyzz, g_y_0_xzzz_xxyyzzz, g_y_0_xzzz_xxyzzz, g_y_0_xzzz_xxyzzzz, g_y_0_xzzz_xxzzzz, g_y_0_xzzz_xxzzzzz, g_y_0_xzzz_xyyyyy, g_y_0_xzzz_xyyyyyy, g_y_0_xzzz_xyyyyyz, g_y_0_xzzz_xyyyyz, g_y_0_xzzz_xyyyyzz, g_y_0_xzzz_xyyyzz, g_y_0_xzzz_xyyyzzz, g_y_0_xzzz_xyyzzz, g_y_0_xzzz_xyyzzzz, g_y_0_xzzz_xyzzzz, g_y_0_xzzz_xyzzzzz, g_y_0_xzzz_xzzzzz, g_y_0_xzzz_xzzzzzz, g_y_0_xzzz_yyyyyy, g_y_0_xzzz_yyyyyz, g_y_0_xzzz_yyyyzz, g_y_0_xzzz_yyyzzz, g_y_0_xzzz_yyzzzz, g_y_0_xzzz_yzzzzz, g_y_0_xzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzz_xxxxxx[k] = -g_y_0_xzzz_xxxxxx[k] * cd_x[k] + g_y_0_xzzz_xxxxxxx[k];

                g_y_0_xxzzz_xxxxxy[k] = -g_y_0_xzzz_xxxxxy[k] * cd_x[k] + g_y_0_xzzz_xxxxxxy[k];

                g_y_0_xxzzz_xxxxxz[k] = -g_y_0_xzzz_xxxxxz[k] * cd_x[k] + g_y_0_xzzz_xxxxxxz[k];

                g_y_0_xxzzz_xxxxyy[k] = -g_y_0_xzzz_xxxxyy[k] * cd_x[k] + g_y_0_xzzz_xxxxxyy[k];

                g_y_0_xxzzz_xxxxyz[k] = -g_y_0_xzzz_xxxxyz[k] * cd_x[k] + g_y_0_xzzz_xxxxxyz[k];

                g_y_0_xxzzz_xxxxzz[k] = -g_y_0_xzzz_xxxxzz[k] * cd_x[k] + g_y_0_xzzz_xxxxxzz[k];

                g_y_0_xxzzz_xxxyyy[k] = -g_y_0_xzzz_xxxyyy[k] * cd_x[k] + g_y_0_xzzz_xxxxyyy[k];

                g_y_0_xxzzz_xxxyyz[k] = -g_y_0_xzzz_xxxyyz[k] * cd_x[k] + g_y_0_xzzz_xxxxyyz[k];

                g_y_0_xxzzz_xxxyzz[k] = -g_y_0_xzzz_xxxyzz[k] * cd_x[k] + g_y_0_xzzz_xxxxyzz[k];

                g_y_0_xxzzz_xxxzzz[k] = -g_y_0_xzzz_xxxzzz[k] * cd_x[k] + g_y_0_xzzz_xxxxzzz[k];

                g_y_0_xxzzz_xxyyyy[k] = -g_y_0_xzzz_xxyyyy[k] * cd_x[k] + g_y_0_xzzz_xxxyyyy[k];

                g_y_0_xxzzz_xxyyyz[k] = -g_y_0_xzzz_xxyyyz[k] * cd_x[k] + g_y_0_xzzz_xxxyyyz[k];

                g_y_0_xxzzz_xxyyzz[k] = -g_y_0_xzzz_xxyyzz[k] * cd_x[k] + g_y_0_xzzz_xxxyyzz[k];

                g_y_0_xxzzz_xxyzzz[k] = -g_y_0_xzzz_xxyzzz[k] * cd_x[k] + g_y_0_xzzz_xxxyzzz[k];

                g_y_0_xxzzz_xxzzzz[k] = -g_y_0_xzzz_xxzzzz[k] * cd_x[k] + g_y_0_xzzz_xxxzzzz[k];

                g_y_0_xxzzz_xyyyyy[k] = -g_y_0_xzzz_xyyyyy[k] * cd_x[k] + g_y_0_xzzz_xxyyyyy[k];

                g_y_0_xxzzz_xyyyyz[k] = -g_y_0_xzzz_xyyyyz[k] * cd_x[k] + g_y_0_xzzz_xxyyyyz[k];

                g_y_0_xxzzz_xyyyzz[k] = -g_y_0_xzzz_xyyyzz[k] * cd_x[k] + g_y_0_xzzz_xxyyyzz[k];

                g_y_0_xxzzz_xyyzzz[k] = -g_y_0_xzzz_xyyzzz[k] * cd_x[k] + g_y_0_xzzz_xxyyzzz[k];

                g_y_0_xxzzz_xyzzzz[k] = -g_y_0_xzzz_xyzzzz[k] * cd_x[k] + g_y_0_xzzz_xxyzzzz[k];

                g_y_0_xxzzz_xzzzzz[k] = -g_y_0_xzzz_xzzzzz[k] * cd_x[k] + g_y_0_xzzz_xxzzzzz[k];

                g_y_0_xxzzz_yyyyyy[k] = -g_y_0_xzzz_yyyyyy[k] * cd_x[k] + g_y_0_xzzz_xyyyyyy[k];

                g_y_0_xxzzz_yyyyyz[k] = -g_y_0_xzzz_yyyyyz[k] * cd_x[k] + g_y_0_xzzz_xyyyyyz[k];

                g_y_0_xxzzz_yyyyzz[k] = -g_y_0_xzzz_yyyyzz[k] * cd_x[k] + g_y_0_xzzz_xyyyyzz[k];

                g_y_0_xxzzz_yyyzzz[k] = -g_y_0_xzzz_yyyzzz[k] * cd_x[k] + g_y_0_xzzz_xyyyzzz[k];

                g_y_0_xxzzz_yyzzzz[k] = -g_y_0_xzzz_yyzzzz[k] * cd_x[k] + g_y_0_xzzz_xyyzzzz[k];

                g_y_0_xxzzz_yzzzzz[k] = -g_y_0_xzzz_yzzzzz[k] * cd_x[k] + g_y_0_xzzz_xyzzzzz[k];

                g_y_0_xxzzz_zzzzzz[k] = -g_y_0_xzzz_zzzzzz[k] * cd_x[k] + g_y_0_xzzz_xzzzzzz[k];
            }

            /// Set up 280-308 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 280);

            auto g_y_0_xyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 281);

            auto g_y_0_xyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 282);

            auto g_y_0_xyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 283);

            auto g_y_0_xyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 284);

            auto g_y_0_xyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 285);

            auto g_y_0_xyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 286);

            auto g_y_0_xyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 287);

            auto g_y_0_xyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 288);

            auto g_y_0_xyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 289);

            auto g_y_0_xyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 290);

            auto g_y_0_xyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 291);

            auto g_y_0_xyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 292);

            auto g_y_0_xyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 293);

            auto g_y_0_xyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 294);

            auto g_y_0_xyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 295);

            auto g_y_0_xyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 296);

            auto g_y_0_xyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 297);

            auto g_y_0_xyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 298);

            auto g_y_0_xyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 299);

            auto g_y_0_xyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 300);

            auto g_y_0_xyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 301);

            auto g_y_0_xyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 302);

            auto g_y_0_xyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 303);

            auto g_y_0_xyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 304);

            auto g_y_0_xyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 305);

            auto g_y_0_xyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 306);

            auto g_y_0_xyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 307);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyy_xxxxxx, g_y_0_xyyyy_xxxxxy, g_y_0_xyyyy_xxxxxz, g_y_0_xyyyy_xxxxyy, g_y_0_xyyyy_xxxxyz, g_y_0_xyyyy_xxxxzz, g_y_0_xyyyy_xxxyyy, g_y_0_xyyyy_xxxyyz, g_y_0_xyyyy_xxxyzz, g_y_0_xyyyy_xxxzzz, g_y_0_xyyyy_xxyyyy, g_y_0_xyyyy_xxyyyz, g_y_0_xyyyy_xxyyzz, g_y_0_xyyyy_xxyzzz, g_y_0_xyyyy_xxzzzz, g_y_0_xyyyy_xyyyyy, g_y_0_xyyyy_xyyyyz, g_y_0_xyyyy_xyyyzz, g_y_0_xyyyy_xyyzzz, g_y_0_xyyyy_xyzzzz, g_y_0_xyyyy_xzzzzz, g_y_0_xyyyy_yyyyyy, g_y_0_xyyyy_yyyyyz, g_y_0_xyyyy_yyyyzz, g_y_0_xyyyy_yyyzzz, g_y_0_xyyyy_yyzzzz, g_y_0_xyyyy_yzzzzz, g_y_0_xyyyy_zzzzzz, g_y_0_yyyy_xxxxxx, g_y_0_yyyy_xxxxxxx, g_y_0_yyyy_xxxxxxy, g_y_0_yyyy_xxxxxxz, g_y_0_yyyy_xxxxxy, g_y_0_yyyy_xxxxxyy, g_y_0_yyyy_xxxxxyz, g_y_0_yyyy_xxxxxz, g_y_0_yyyy_xxxxxzz, g_y_0_yyyy_xxxxyy, g_y_0_yyyy_xxxxyyy, g_y_0_yyyy_xxxxyyz, g_y_0_yyyy_xxxxyz, g_y_0_yyyy_xxxxyzz, g_y_0_yyyy_xxxxzz, g_y_0_yyyy_xxxxzzz, g_y_0_yyyy_xxxyyy, g_y_0_yyyy_xxxyyyy, g_y_0_yyyy_xxxyyyz, g_y_0_yyyy_xxxyyz, g_y_0_yyyy_xxxyyzz, g_y_0_yyyy_xxxyzz, g_y_0_yyyy_xxxyzzz, g_y_0_yyyy_xxxzzz, g_y_0_yyyy_xxxzzzz, g_y_0_yyyy_xxyyyy, g_y_0_yyyy_xxyyyyy, g_y_0_yyyy_xxyyyyz, g_y_0_yyyy_xxyyyz, g_y_0_yyyy_xxyyyzz, g_y_0_yyyy_xxyyzz, g_y_0_yyyy_xxyyzzz, g_y_0_yyyy_xxyzzz, g_y_0_yyyy_xxyzzzz, g_y_0_yyyy_xxzzzz, g_y_0_yyyy_xxzzzzz, g_y_0_yyyy_xyyyyy, g_y_0_yyyy_xyyyyyy, g_y_0_yyyy_xyyyyyz, g_y_0_yyyy_xyyyyz, g_y_0_yyyy_xyyyyzz, g_y_0_yyyy_xyyyzz, g_y_0_yyyy_xyyyzzz, g_y_0_yyyy_xyyzzz, g_y_0_yyyy_xyyzzzz, g_y_0_yyyy_xyzzzz, g_y_0_yyyy_xyzzzzz, g_y_0_yyyy_xzzzzz, g_y_0_yyyy_xzzzzzz, g_y_0_yyyy_yyyyyy, g_y_0_yyyy_yyyyyz, g_y_0_yyyy_yyyyzz, g_y_0_yyyy_yyyzzz, g_y_0_yyyy_yyzzzz, g_y_0_yyyy_yzzzzz, g_y_0_yyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyy_xxxxxx[k] = -g_y_0_yyyy_xxxxxx[k] * cd_x[k] + g_y_0_yyyy_xxxxxxx[k];

                g_y_0_xyyyy_xxxxxy[k] = -g_y_0_yyyy_xxxxxy[k] * cd_x[k] + g_y_0_yyyy_xxxxxxy[k];

                g_y_0_xyyyy_xxxxxz[k] = -g_y_0_yyyy_xxxxxz[k] * cd_x[k] + g_y_0_yyyy_xxxxxxz[k];

                g_y_0_xyyyy_xxxxyy[k] = -g_y_0_yyyy_xxxxyy[k] * cd_x[k] + g_y_0_yyyy_xxxxxyy[k];

                g_y_0_xyyyy_xxxxyz[k] = -g_y_0_yyyy_xxxxyz[k] * cd_x[k] + g_y_0_yyyy_xxxxxyz[k];

                g_y_0_xyyyy_xxxxzz[k] = -g_y_0_yyyy_xxxxzz[k] * cd_x[k] + g_y_0_yyyy_xxxxxzz[k];

                g_y_0_xyyyy_xxxyyy[k] = -g_y_0_yyyy_xxxyyy[k] * cd_x[k] + g_y_0_yyyy_xxxxyyy[k];

                g_y_0_xyyyy_xxxyyz[k] = -g_y_0_yyyy_xxxyyz[k] * cd_x[k] + g_y_0_yyyy_xxxxyyz[k];

                g_y_0_xyyyy_xxxyzz[k] = -g_y_0_yyyy_xxxyzz[k] * cd_x[k] + g_y_0_yyyy_xxxxyzz[k];

                g_y_0_xyyyy_xxxzzz[k] = -g_y_0_yyyy_xxxzzz[k] * cd_x[k] + g_y_0_yyyy_xxxxzzz[k];

                g_y_0_xyyyy_xxyyyy[k] = -g_y_0_yyyy_xxyyyy[k] * cd_x[k] + g_y_0_yyyy_xxxyyyy[k];

                g_y_0_xyyyy_xxyyyz[k] = -g_y_0_yyyy_xxyyyz[k] * cd_x[k] + g_y_0_yyyy_xxxyyyz[k];

                g_y_0_xyyyy_xxyyzz[k] = -g_y_0_yyyy_xxyyzz[k] * cd_x[k] + g_y_0_yyyy_xxxyyzz[k];

                g_y_0_xyyyy_xxyzzz[k] = -g_y_0_yyyy_xxyzzz[k] * cd_x[k] + g_y_0_yyyy_xxxyzzz[k];

                g_y_0_xyyyy_xxzzzz[k] = -g_y_0_yyyy_xxzzzz[k] * cd_x[k] + g_y_0_yyyy_xxxzzzz[k];

                g_y_0_xyyyy_xyyyyy[k] = -g_y_0_yyyy_xyyyyy[k] * cd_x[k] + g_y_0_yyyy_xxyyyyy[k];

                g_y_0_xyyyy_xyyyyz[k] = -g_y_0_yyyy_xyyyyz[k] * cd_x[k] + g_y_0_yyyy_xxyyyyz[k];

                g_y_0_xyyyy_xyyyzz[k] = -g_y_0_yyyy_xyyyzz[k] * cd_x[k] + g_y_0_yyyy_xxyyyzz[k];

                g_y_0_xyyyy_xyyzzz[k] = -g_y_0_yyyy_xyyzzz[k] * cd_x[k] + g_y_0_yyyy_xxyyzzz[k];

                g_y_0_xyyyy_xyzzzz[k] = -g_y_0_yyyy_xyzzzz[k] * cd_x[k] + g_y_0_yyyy_xxyzzzz[k];

                g_y_0_xyyyy_xzzzzz[k] = -g_y_0_yyyy_xzzzzz[k] * cd_x[k] + g_y_0_yyyy_xxzzzzz[k];

                g_y_0_xyyyy_yyyyyy[k] = -g_y_0_yyyy_yyyyyy[k] * cd_x[k] + g_y_0_yyyy_xyyyyyy[k];

                g_y_0_xyyyy_yyyyyz[k] = -g_y_0_yyyy_yyyyyz[k] * cd_x[k] + g_y_0_yyyy_xyyyyyz[k];

                g_y_0_xyyyy_yyyyzz[k] = -g_y_0_yyyy_yyyyzz[k] * cd_x[k] + g_y_0_yyyy_xyyyyzz[k];

                g_y_0_xyyyy_yyyzzz[k] = -g_y_0_yyyy_yyyzzz[k] * cd_x[k] + g_y_0_yyyy_xyyyzzz[k];

                g_y_0_xyyyy_yyzzzz[k] = -g_y_0_yyyy_yyzzzz[k] * cd_x[k] + g_y_0_yyyy_xyyzzzz[k];

                g_y_0_xyyyy_yzzzzz[k] = -g_y_0_yyyy_yzzzzz[k] * cd_x[k] + g_y_0_yyyy_xyzzzzz[k];

                g_y_0_xyyyy_zzzzzz[k] = -g_y_0_yyyy_zzzzzz[k] * cd_x[k] + g_y_0_yyyy_xzzzzzz[k];
            }

            /// Set up 308-336 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 308);

            auto g_y_0_xyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 309);

            auto g_y_0_xyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 310);

            auto g_y_0_xyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 311);

            auto g_y_0_xyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 312);

            auto g_y_0_xyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 313);

            auto g_y_0_xyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 314);

            auto g_y_0_xyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 315);

            auto g_y_0_xyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 316);

            auto g_y_0_xyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 317);

            auto g_y_0_xyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 318);

            auto g_y_0_xyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 319);

            auto g_y_0_xyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 320);

            auto g_y_0_xyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 321);

            auto g_y_0_xyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 322);

            auto g_y_0_xyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 323);

            auto g_y_0_xyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 324);

            auto g_y_0_xyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 325);

            auto g_y_0_xyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 326);

            auto g_y_0_xyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 327);

            auto g_y_0_xyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 328);

            auto g_y_0_xyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 329);

            auto g_y_0_xyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 330);

            auto g_y_0_xyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 331);

            auto g_y_0_xyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 332);

            auto g_y_0_xyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 333);

            auto g_y_0_xyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 334);

            auto g_y_0_xyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 335);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyz_xxxxxx, g_y_0_xyyyz_xxxxxy, g_y_0_xyyyz_xxxxxz, g_y_0_xyyyz_xxxxyy, g_y_0_xyyyz_xxxxyz, g_y_0_xyyyz_xxxxzz, g_y_0_xyyyz_xxxyyy, g_y_0_xyyyz_xxxyyz, g_y_0_xyyyz_xxxyzz, g_y_0_xyyyz_xxxzzz, g_y_0_xyyyz_xxyyyy, g_y_0_xyyyz_xxyyyz, g_y_0_xyyyz_xxyyzz, g_y_0_xyyyz_xxyzzz, g_y_0_xyyyz_xxzzzz, g_y_0_xyyyz_xyyyyy, g_y_0_xyyyz_xyyyyz, g_y_0_xyyyz_xyyyzz, g_y_0_xyyyz_xyyzzz, g_y_0_xyyyz_xyzzzz, g_y_0_xyyyz_xzzzzz, g_y_0_xyyyz_yyyyyy, g_y_0_xyyyz_yyyyyz, g_y_0_xyyyz_yyyyzz, g_y_0_xyyyz_yyyzzz, g_y_0_xyyyz_yyzzzz, g_y_0_xyyyz_yzzzzz, g_y_0_xyyyz_zzzzzz, g_y_0_yyyz_xxxxxx, g_y_0_yyyz_xxxxxxx, g_y_0_yyyz_xxxxxxy, g_y_0_yyyz_xxxxxxz, g_y_0_yyyz_xxxxxy, g_y_0_yyyz_xxxxxyy, g_y_0_yyyz_xxxxxyz, g_y_0_yyyz_xxxxxz, g_y_0_yyyz_xxxxxzz, g_y_0_yyyz_xxxxyy, g_y_0_yyyz_xxxxyyy, g_y_0_yyyz_xxxxyyz, g_y_0_yyyz_xxxxyz, g_y_0_yyyz_xxxxyzz, g_y_0_yyyz_xxxxzz, g_y_0_yyyz_xxxxzzz, g_y_0_yyyz_xxxyyy, g_y_0_yyyz_xxxyyyy, g_y_0_yyyz_xxxyyyz, g_y_0_yyyz_xxxyyz, g_y_0_yyyz_xxxyyzz, g_y_0_yyyz_xxxyzz, g_y_0_yyyz_xxxyzzz, g_y_0_yyyz_xxxzzz, g_y_0_yyyz_xxxzzzz, g_y_0_yyyz_xxyyyy, g_y_0_yyyz_xxyyyyy, g_y_0_yyyz_xxyyyyz, g_y_0_yyyz_xxyyyz, g_y_0_yyyz_xxyyyzz, g_y_0_yyyz_xxyyzz, g_y_0_yyyz_xxyyzzz, g_y_0_yyyz_xxyzzz, g_y_0_yyyz_xxyzzzz, g_y_0_yyyz_xxzzzz, g_y_0_yyyz_xxzzzzz, g_y_0_yyyz_xyyyyy, g_y_0_yyyz_xyyyyyy, g_y_0_yyyz_xyyyyyz, g_y_0_yyyz_xyyyyz, g_y_0_yyyz_xyyyyzz, g_y_0_yyyz_xyyyzz, g_y_0_yyyz_xyyyzzz, g_y_0_yyyz_xyyzzz, g_y_0_yyyz_xyyzzzz, g_y_0_yyyz_xyzzzz, g_y_0_yyyz_xyzzzzz, g_y_0_yyyz_xzzzzz, g_y_0_yyyz_xzzzzzz, g_y_0_yyyz_yyyyyy, g_y_0_yyyz_yyyyyz, g_y_0_yyyz_yyyyzz, g_y_0_yyyz_yyyzzz, g_y_0_yyyz_yyzzzz, g_y_0_yyyz_yzzzzz, g_y_0_yyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyz_xxxxxx[k] = -g_y_0_yyyz_xxxxxx[k] * cd_x[k] + g_y_0_yyyz_xxxxxxx[k];

                g_y_0_xyyyz_xxxxxy[k] = -g_y_0_yyyz_xxxxxy[k] * cd_x[k] + g_y_0_yyyz_xxxxxxy[k];

                g_y_0_xyyyz_xxxxxz[k] = -g_y_0_yyyz_xxxxxz[k] * cd_x[k] + g_y_0_yyyz_xxxxxxz[k];

                g_y_0_xyyyz_xxxxyy[k] = -g_y_0_yyyz_xxxxyy[k] * cd_x[k] + g_y_0_yyyz_xxxxxyy[k];

                g_y_0_xyyyz_xxxxyz[k] = -g_y_0_yyyz_xxxxyz[k] * cd_x[k] + g_y_0_yyyz_xxxxxyz[k];

                g_y_0_xyyyz_xxxxzz[k] = -g_y_0_yyyz_xxxxzz[k] * cd_x[k] + g_y_0_yyyz_xxxxxzz[k];

                g_y_0_xyyyz_xxxyyy[k] = -g_y_0_yyyz_xxxyyy[k] * cd_x[k] + g_y_0_yyyz_xxxxyyy[k];

                g_y_0_xyyyz_xxxyyz[k] = -g_y_0_yyyz_xxxyyz[k] * cd_x[k] + g_y_0_yyyz_xxxxyyz[k];

                g_y_0_xyyyz_xxxyzz[k] = -g_y_0_yyyz_xxxyzz[k] * cd_x[k] + g_y_0_yyyz_xxxxyzz[k];

                g_y_0_xyyyz_xxxzzz[k] = -g_y_0_yyyz_xxxzzz[k] * cd_x[k] + g_y_0_yyyz_xxxxzzz[k];

                g_y_0_xyyyz_xxyyyy[k] = -g_y_0_yyyz_xxyyyy[k] * cd_x[k] + g_y_0_yyyz_xxxyyyy[k];

                g_y_0_xyyyz_xxyyyz[k] = -g_y_0_yyyz_xxyyyz[k] * cd_x[k] + g_y_0_yyyz_xxxyyyz[k];

                g_y_0_xyyyz_xxyyzz[k] = -g_y_0_yyyz_xxyyzz[k] * cd_x[k] + g_y_0_yyyz_xxxyyzz[k];

                g_y_0_xyyyz_xxyzzz[k] = -g_y_0_yyyz_xxyzzz[k] * cd_x[k] + g_y_0_yyyz_xxxyzzz[k];

                g_y_0_xyyyz_xxzzzz[k] = -g_y_0_yyyz_xxzzzz[k] * cd_x[k] + g_y_0_yyyz_xxxzzzz[k];

                g_y_0_xyyyz_xyyyyy[k] = -g_y_0_yyyz_xyyyyy[k] * cd_x[k] + g_y_0_yyyz_xxyyyyy[k];

                g_y_0_xyyyz_xyyyyz[k] = -g_y_0_yyyz_xyyyyz[k] * cd_x[k] + g_y_0_yyyz_xxyyyyz[k];

                g_y_0_xyyyz_xyyyzz[k] = -g_y_0_yyyz_xyyyzz[k] * cd_x[k] + g_y_0_yyyz_xxyyyzz[k];

                g_y_0_xyyyz_xyyzzz[k] = -g_y_0_yyyz_xyyzzz[k] * cd_x[k] + g_y_0_yyyz_xxyyzzz[k];

                g_y_0_xyyyz_xyzzzz[k] = -g_y_0_yyyz_xyzzzz[k] * cd_x[k] + g_y_0_yyyz_xxyzzzz[k];

                g_y_0_xyyyz_xzzzzz[k] = -g_y_0_yyyz_xzzzzz[k] * cd_x[k] + g_y_0_yyyz_xxzzzzz[k];

                g_y_0_xyyyz_yyyyyy[k] = -g_y_0_yyyz_yyyyyy[k] * cd_x[k] + g_y_0_yyyz_xyyyyyy[k];

                g_y_0_xyyyz_yyyyyz[k] = -g_y_0_yyyz_yyyyyz[k] * cd_x[k] + g_y_0_yyyz_xyyyyyz[k];

                g_y_0_xyyyz_yyyyzz[k] = -g_y_0_yyyz_yyyyzz[k] * cd_x[k] + g_y_0_yyyz_xyyyyzz[k];

                g_y_0_xyyyz_yyyzzz[k] = -g_y_0_yyyz_yyyzzz[k] * cd_x[k] + g_y_0_yyyz_xyyyzzz[k];

                g_y_0_xyyyz_yyzzzz[k] = -g_y_0_yyyz_yyzzzz[k] * cd_x[k] + g_y_0_yyyz_xyyzzzz[k];

                g_y_0_xyyyz_yzzzzz[k] = -g_y_0_yyyz_yzzzzz[k] * cd_x[k] + g_y_0_yyyz_xyzzzzz[k];

                g_y_0_xyyyz_zzzzzz[k] = -g_y_0_yyyz_zzzzzz[k] * cd_x[k] + g_y_0_yyyz_xzzzzzz[k];
            }

            /// Set up 336-364 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 336);

            auto g_y_0_xyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 337);

            auto g_y_0_xyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 338);

            auto g_y_0_xyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 339);

            auto g_y_0_xyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 340);

            auto g_y_0_xyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 341);

            auto g_y_0_xyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 342);

            auto g_y_0_xyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 343);

            auto g_y_0_xyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 344);

            auto g_y_0_xyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 345);

            auto g_y_0_xyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 346);

            auto g_y_0_xyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 347);

            auto g_y_0_xyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 348);

            auto g_y_0_xyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 349);

            auto g_y_0_xyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 350);

            auto g_y_0_xyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 351);

            auto g_y_0_xyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 352);

            auto g_y_0_xyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 353);

            auto g_y_0_xyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 354);

            auto g_y_0_xyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 355);

            auto g_y_0_xyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 356);

            auto g_y_0_xyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 357);

            auto g_y_0_xyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 358);

            auto g_y_0_xyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 359);

            auto g_y_0_xyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 360);

            auto g_y_0_xyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 361);

            auto g_y_0_xyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 362);

            auto g_y_0_xyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 363);

            #pragma omp simd aligned(cd_x, g_y_0_xyyzz_xxxxxx, g_y_0_xyyzz_xxxxxy, g_y_0_xyyzz_xxxxxz, g_y_0_xyyzz_xxxxyy, g_y_0_xyyzz_xxxxyz, g_y_0_xyyzz_xxxxzz, g_y_0_xyyzz_xxxyyy, g_y_0_xyyzz_xxxyyz, g_y_0_xyyzz_xxxyzz, g_y_0_xyyzz_xxxzzz, g_y_0_xyyzz_xxyyyy, g_y_0_xyyzz_xxyyyz, g_y_0_xyyzz_xxyyzz, g_y_0_xyyzz_xxyzzz, g_y_0_xyyzz_xxzzzz, g_y_0_xyyzz_xyyyyy, g_y_0_xyyzz_xyyyyz, g_y_0_xyyzz_xyyyzz, g_y_0_xyyzz_xyyzzz, g_y_0_xyyzz_xyzzzz, g_y_0_xyyzz_xzzzzz, g_y_0_xyyzz_yyyyyy, g_y_0_xyyzz_yyyyyz, g_y_0_xyyzz_yyyyzz, g_y_0_xyyzz_yyyzzz, g_y_0_xyyzz_yyzzzz, g_y_0_xyyzz_yzzzzz, g_y_0_xyyzz_zzzzzz, g_y_0_yyzz_xxxxxx, g_y_0_yyzz_xxxxxxx, g_y_0_yyzz_xxxxxxy, g_y_0_yyzz_xxxxxxz, g_y_0_yyzz_xxxxxy, g_y_0_yyzz_xxxxxyy, g_y_0_yyzz_xxxxxyz, g_y_0_yyzz_xxxxxz, g_y_0_yyzz_xxxxxzz, g_y_0_yyzz_xxxxyy, g_y_0_yyzz_xxxxyyy, g_y_0_yyzz_xxxxyyz, g_y_0_yyzz_xxxxyz, g_y_0_yyzz_xxxxyzz, g_y_0_yyzz_xxxxzz, g_y_0_yyzz_xxxxzzz, g_y_0_yyzz_xxxyyy, g_y_0_yyzz_xxxyyyy, g_y_0_yyzz_xxxyyyz, g_y_0_yyzz_xxxyyz, g_y_0_yyzz_xxxyyzz, g_y_0_yyzz_xxxyzz, g_y_0_yyzz_xxxyzzz, g_y_0_yyzz_xxxzzz, g_y_0_yyzz_xxxzzzz, g_y_0_yyzz_xxyyyy, g_y_0_yyzz_xxyyyyy, g_y_0_yyzz_xxyyyyz, g_y_0_yyzz_xxyyyz, g_y_0_yyzz_xxyyyzz, g_y_0_yyzz_xxyyzz, g_y_0_yyzz_xxyyzzz, g_y_0_yyzz_xxyzzz, g_y_0_yyzz_xxyzzzz, g_y_0_yyzz_xxzzzz, g_y_0_yyzz_xxzzzzz, g_y_0_yyzz_xyyyyy, g_y_0_yyzz_xyyyyyy, g_y_0_yyzz_xyyyyyz, g_y_0_yyzz_xyyyyz, g_y_0_yyzz_xyyyyzz, g_y_0_yyzz_xyyyzz, g_y_0_yyzz_xyyyzzz, g_y_0_yyzz_xyyzzz, g_y_0_yyzz_xyyzzzz, g_y_0_yyzz_xyzzzz, g_y_0_yyzz_xyzzzzz, g_y_0_yyzz_xzzzzz, g_y_0_yyzz_xzzzzzz, g_y_0_yyzz_yyyyyy, g_y_0_yyzz_yyyyyz, g_y_0_yyzz_yyyyzz, g_y_0_yyzz_yyyzzz, g_y_0_yyzz_yyzzzz, g_y_0_yyzz_yzzzzz, g_y_0_yyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzz_xxxxxx[k] = -g_y_0_yyzz_xxxxxx[k] * cd_x[k] + g_y_0_yyzz_xxxxxxx[k];

                g_y_0_xyyzz_xxxxxy[k] = -g_y_0_yyzz_xxxxxy[k] * cd_x[k] + g_y_0_yyzz_xxxxxxy[k];

                g_y_0_xyyzz_xxxxxz[k] = -g_y_0_yyzz_xxxxxz[k] * cd_x[k] + g_y_0_yyzz_xxxxxxz[k];

                g_y_0_xyyzz_xxxxyy[k] = -g_y_0_yyzz_xxxxyy[k] * cd_x[k] + g_y_0_yyzz_xxxxxyy[k];

                g_y_0_xyyzz_xxxxyz[k] = -g_y_0_yyzz_xxxxyz[k] * cd_x[k] + g_y_0_yyzz_xxxxxyz[k];

                g_y_0_xyyzz_xxxxzz[k] = -g_y_0_yyzz_xxxxzz[k] * cd_x[k] + g_y_0_yyzz_xxxxxzz[k];

                g_y_0_xyyzz_xxxyyy[k] = -g_y_0_yyzz_xxxyyy[k] * cd_x[k] + g_y_0_yyzz_xxxxyyy[k];

                g_y_0_xyyzz_xxxyyz[k] = -g_y_0_yyzz_xxxyyz[k] * cd_x[k] + g_y_0_yyzz_xxxxyyz[k];

                g_y_0_xyyzz_xxxyzz[k] = -g_y_0_yyzz_xxxyzz[k] * cd_x[k] + g_y_0_yyzz_xxxxyzz[k];

                g_y_0_xyyzz_xxxzzz[k] = -g_y_0_yyzz_xxxzzz[k] * cd_x[k] + g_y_0_yyzz_xxxxzzz[k];

                g_y_0_xyyzz_xxyyyy[k] = -g_y_0_yyzz_xxyyyy[k] * cd_x[k] + g_y_0_yyzz_xxxyyyy[k];

                g_y_0_xyyzz_xxyyyz[k] = -g_y_0_yyzz_xxyyyz[k] * cd_x[k] + g_y_0_yyzz_xxxyyyz[k];

                g_y_0_xyyzz_xxyyzz[k] = -g_y_0_yyzz_xxyyzz[k] * cd_x[k] + g_y_0_yyzz_xxxyyzz[k];

                g_y_0_xyyzz_xxyzzz[k] = -g_y_0_yyzz_xxyzzz[k] * cd_x[k] + g_y_0_yyzz_xxxyzzz[k];

                g_y_0_xyyzz_xxzzzz[k] = -g_y_0_yyzz_xxzzzz[k] * cd_x[k] + g_y_0_yyzz_xxxzzzz[k];

                g_y_0_xyyzz_xyyyyy[k] = -g_y_0_yyzz_xyyyyy[k] * cd_x[k] + g_y_0_yyzz_xxyyyyy[k];

                g_y_0_xyyzz_xyyyyz[k] = -g_y_0_yyzz_xyyyyz[k] * cd_x[k] + g_y_0_yyzz_xxyyyyz[k];

                g_y_0_xyyzz_xyyyzz[k] = -g_y_0_yyzz_xyyyzz[k] * cd_x[k] + g_y_0_yyzz_xxyyyzz[k];

                g_y_0_xyyzz_xyyzzz[k] = -g_y_0_yyzz_xyyzzz[k] * cd_x[k] + g_y_0_yyzz_xxyyzzz[k];

                g_y_0_xyyzz_xyzzzz[k] = -g_y_0_yyzz_xyzzzz[k] * cd_x[k] + g_y_0_yyzz_xxyzzzz[k];

                g_y_0_xyyzz_xzzzzz[k] = -g_y_0_yyzz_xzzzzz[k] * cd_x[k] + g_y_0_yyzz_xxzzzzz[k];

                g_y_0_xyyzz_yyyyyy[k] = -g_y_0_yyzz_yyyyyy[k] * cd_x[k] + g_y_0_yyzz_xyyyyyy[k];

                g_y_0_xyyzz_yyyyyz[k] = -g_y_0_yyzz_yyyyyz[k] * cd_x[k] + g_y_0_yyzz_xyyyyyz[k];

                g_y_0_xyyzz_yyyyzz[k] = -g_y_0_yyzz_yyyyzz[k] * cd_x[k] + g_y_0_yyzz_xyyyyzz[k];

                g_y_0_xyyzz_yyyzzz[k] = -g_y_0_yyzz_yyyzzz[k] * cd_x[k] + g_y_0_yyzz_xyyyzzz[k];

                g_y_0_xyyzz_yyzzzz[k] = -g_y_0_yyzz_yyzzzz[k] * cd_x[k] + g_y_0_yyzz_xyyzzzz[k];

                g_y_0_xyyzz_yzzzzz[k] = -g_y_0_yyzz_yzzzzz[k] * cd_x[k] + g_y_0_yyzz_xyzzzzz[k];

                g_y_0_xyyzz_zzzzzz[k] = -g_y_0_yyzz_zzzzzz[k] * cd_x[k] + g_y_0_yyzz_xzzzzzz[k];
            }

            /// Set up 364-392 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 364);

            auto g_y_0_xyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 365);

            auto g_y_0_xyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 366);

            auto g_y_0_xyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 367);

            auto g_y_0_xyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 368);

            auto g_y_0_xyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 369);

            auto g_y_0_xyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 370);

            auto g_y_0_xyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 371);

            auto g_y_0_xyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 372);

            auto g_y_0_xyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 373);

            auto g_y_0_xyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 374);

            auto g_y_0_xyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 375);

            auto g_y_0_xyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 376);

            auto g_y_0_xyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 377);

            auto g_y_0_xyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 378);

            auto g_y_0_xyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 379);

            auto g_y_0_xyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 380);

            auto g_y_0_xyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 381);

            auto g_y_0_xyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 382);

            auto g_y_0_xyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 383);

            auto g_y_0_xyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 384);

            auto g_y_0_xyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 385);

            auto g_y_0_xyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 386);

            auto g_y_0_xyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 387);

            auto g_y_0_xyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 388);

            auto g_y_0_xyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 389);

            auto g_y_0_xyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 390);

            auto g_y_0_xyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 391);

            #pragma omp simd aligned(cd_x, g_y_0_xyzzz_xxxxxx, g_y_0_xyzzz_xxxxxy, g_y_0_xyzzz_xxxxxz, g_y_0_xyzzz_xxxxyy, g_y_0_xyzzz_xxxxyz, g_y_0_xyzzz_xxxxzz, g_y_0_xyzzz_xxxyyy, g_y_0_xyzzz_xxxyyz, g_y_0_xyzzz_xxxyzz, g_y_0_xyzzz_xxxzzz, g_y_0_xyzzz_xxyyyy, g_y_0_xyzzz_xxyyyz, g_y_0_xyzzz_xxyyzz, g_y_0_xyzzz_xxyzzz, g_y_0_xyzzz_xxzzzz, g_y_0_xyzzz_xyyyyy, g_y_0_xyzzz_xyyyyz, g_y_0_xyzzz_xyyyzz, g_y_0_xyzzz_xyyzzz, g_y_0_xyzzz_xyzzzz, g_y_0_xyzzz_xzzzzz, g_y_0_xyzzz_yyyyyy, g_y_0_xyzzz_yyyyyz, g_y_0_xyzzz_yyyyzz, g_y_0_xyzzz_yyyzzz, g_y_0_xyzzz_yyzzzz, g_y_0_xyzzz_yzzzzz, g_y_0_xyzzz_zzzzzz, g_y_0_yzzz_xxxxxx, g_y_0_yzzz_xxxxxxx, g_y_0_yzzz_xxxxxxy, g_y_0_yzzz_xxxxxxz, g_y_0_yzzz_xxxxxy, g_y_0_yzzz_xxxxxyy, g_y_0_yzzz_xxxxxyz, g_y_0_yzzz_xxxxxz, g_y_0_yzzz_xxxxxzz, g_y_0_yzzz_xxxxyy, g_y_0_yzzz_xxxxyyy, g_y_0_yzzz_xxxxyyz, g_y_0_yzzz_xxxxyz, g_y_0_yzzz_xxxxyzz, g_y_0_yzzz_xxxxzz, g_y_0_yzzz_xxxxzzz, g_y_0_yzzz_xxxyyy, g_y_0_yzzz_xxxyyyy, g_y_0_yzzz_xxxyyyz, g_y_0_yzzz_xxxyyz, g_y_0_yzzz_xxxyyzz, g_y_0_yzzz_xxxyzz, g_y_0_yzzz_xxxyzzz, g_y_0_yzzz_xxxzzz, g_y_0_yzzz_xxxzzzz, g_y_0_yzzz_xxyyyy, g_y_0_yzzz_xxyyyyy, g_y_0_yzzz_xxyyyyz, g_y_0_yzzz_xxyyyz, g_y_0_yzzz_xxyyyzz, g_y_0_yzzz_xxyyzz, g_y_0_yzzz_xxyyzzz, g_y_0_yzzz_xxyzzz, g_y_0_yzzz_xxyzzzz, g_y_0_yzzz_xxzzzz, g_y_0_yzzz_xxzzzzz, g_y_0_yzzz_xyyyyy, g_y_0_yzzz_xyyyyyy, g_y_0_yzzz_xyyyyyz, g_y_0_yzzz_xyyyyz, g_y_0_yzzz_xyyyyzz, g_y_0_yzzz_xyyyzz, g_y_0_yzzz_xyyyzzz, g_y_0_yzzz_xyyzzz, g_y_0_yzzz_xyyzzzz, g_y_0_yzzz_xyzzzz, g_y_0_yzzz_xyzzzzz, g_y_0_yzzz_xzzzzz, g_y_0_yzzz_xzzzzzz, g_y_0_yzzz_yyyyyy, g_y_0_yzzz_yyyyyz, g_y_0_yzzz_yyyyzz, g_y_0_yzzz_yyyzzz, g_y_0_yzzz_yyzzzz, g_y_0_yzzz_yzzzzz, g_y_0_yzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzz_xxxxxx[k] = -g_y_0_yzzz_xxxxxx[k] * cd_x[k] + g_y_0_yzzz_xxxxxxx[k];

                g_y_0_xyzzz_xxxxxy[k] = -g_y_0_yzzz_xxxxxy[k] * cd_x[k] + g_y_0_yzzz_xxxxxxy[k];

                g_y_0_xyzzz_xxxxxz[k] = -g_y_0_yzzz_xxxxxz[k] * cd_x[k] + g_y_0_yzzz_xxxxxxz[k];

                g_y_0_xyzzz_xxxxyy[k] = -g_y_0_yzzz_xxxxyy[k] * cd_x[k] + g_y_0_yzzz_xxxxxyy[k];

                g_y_0_xyzzz_xxxxyz[k] = -g_y_0_yzzz_xxxxyz[k] * cd_x[k] + g_y_0_yzzz_xxxxxyz[k];

                g_y_0_xyzzz_xxxxzz[k] = -g_y_0_yzzz_xxxxzz[k] * cd_x[k] + g_y_0_yzzz_xxxxxzz[k];

                g_y_0_xyzzz_xxxyyy[k] = -g_y_0_yzzz_xxxyyy[k] * cd_x[k] + g_y_0_yzzz_xxxxyyy[k];

                g_y_0_xyzzz_xxxyyz[k] = -g_y_0_yzzz_xxxyyz[k] * cd_x[k] + g_y_0_yzzz_xxxxyyz[k];

                g_y_0_xyzzz_xxxyzz[k] = -g_y_0_yzzz_xxxyzz[k] * cd_x[k] + g_y_0_yzzz_xxxxyzz[k];

                g_y_0_xyzzz_xxxzzz[k] = -g_y_0_yzzz_xxxzzz[k] * cd_x[k] + g_y_0_yzzz_xxxxzzz[k];

                g_y_0_xyzzz_xxyyyy[k] = -g_y_0_yzzz_xxyyyy[k] * cd_x[k] + g_y_0_yzzz_xxxyyyy[k];

                g_y_0_xyzzz_xxyyyz[k] = -g_y_0_yzzz_xxyyyz[k] * cd_x[k] + g_y_0_yzzz_xxxyyyz[k];

                g_y_0_xyzzz_xxyyzz[k] = -g_y_0_yzzz_xxyyzz[k] * cd_x[k] + g_y_0_yzzz_xxxyyzz[k];

                g_y_0_xyzzz_xxyzzz[k] = -g_y_0_yzzz_xxyzzz[k] * cd_x[k] + g_y_0_yzzz_xxxyzzz[k];

                g_y_0_xyzzz_xxzzzz[k] = -g_y_0_yzzz_xxzzzz[k] * cd_x[k] + g_y_0_yzzz_xxxzzzz[k];

                g_y_0_xyzzz_xyyyyy[k] = -g_y_0_yzzz_xyyyyy[k] * cd_x[k] + g_y_0_yzzz_xxyyyyy[k];

                g_y_0_xyzzz_xyyyyz[k] = -g_y_0_yzzz_xyyyyz[k] * cd_x[k] + g_y_0_yzzz_xxyyyyz[k];

                g_y_0_xyzzz_xyyyzz[k] = -g_y_0_yzzz_xyyyzz[k] * cd_x[k] + g_y_0_yzzz_xxyyyzz[k];

                g_y_0_xyzzz_xyyzzz[k] = -g_y_0_yzzz_xyyzzz[k] * cd_x[k] + g_y_0_yzzz_xxyyzzz[k];

                g_y_0_xyzzz_xyzzzz[k] = -g_y_0_yzzz_xyzzzz[k] * cd_x[k] + g_y_0_yzzz_xxyzzzz[k];

                g_y_0_xyzzz_xzzzzz[k] = -g_y_0_yzzz_xzzzzz[k] * cd_x[k] + g_y_0_yzzz_xxzzzzz[k];

                g_y_0_xyzzz_yyyyyy[k] = -g_y_0_yzzz_yyyyyy[k] * cd_x[k] + g_y_0_yzzz_xyyyyyy[k];

                g_y_0_xyzzz_yyyyyz[k] = -g_y_0_yzzz_yyyyyz[k] * cd_x[k] + g_y_0_yzzz_xyyyyyz[k];

                g_y_0_xyzzz_yyyyzz[k] = -g_y_0_yzzz_yyyyzz[k] * cd_x[k] + g_y_0_yzzz_xyyyyzz[k];

                g_y_0_xyzzz_yyyzzz[k] = -g_y_0_yzzz_yyyzzz[k] * cd_x[k] + g_y_0_yzzz_xyyyzzz[k];

                g_y_0_xyzzz_yyzzzz[k] = -g_y_0_yzzz_yyzzzz[k] * cd_x[k] + g_y_0_yzzz_xyyzzzz[k];

                g_y_0_xyzzz_yzzzzz[k] = -g_y_0_yzzz_yzzzzz[k] * cd_x[k] + g_y_0_yzzz_xyzzzzz[k];

                g_y_0_xyzzz_zzzzzz[k] = -g_y_0_yzzz_zzzzzz[k] * cd_x[k] + g_y_0_yzzz_xzzzzzz[k];
            }

            /// Set up 392-420 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 392);

            auto g_y_0_xzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 393);

            auto g_y_0_xzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 394);

            auto g_y_0_xzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 395);

            auto g_y_0_xzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 396);

            auto g_y_0_xzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 397);

            auto g_y_0_xzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 398);

            auto g_y_0_xzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 399);

            auto g_y_0_xzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 400);

            auto g_y_0_xzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 401);

            auto g_y_0_xzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 402);

            auto g_y_0_xzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 403);

            auto g_y_0_xzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 404);

            auto g_y_0_xzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 405);

            auto g_y_0_xzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 406);

            auto g_y_0_xzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 407);

            auto g_y_0_xzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 408);

            auto g_y_0_xzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 409);

            auto g_y_0_xzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 410);

            auto g_y_0_xzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 411);

            auto g_y_0_xzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 412);

            auto g_y_0_xzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 413);

            auto g_y_0_xzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 414);

            auto g_y_0_xzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 415);

            auto g_y_0_xzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 416);

            auto g_y_0_xzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 417);

            auto g_y_0_xzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 418);

            auto g_y_0_xzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 419);

            #pragma omp simd aligned(cd_x, g_y_0_xzzzz_xxxxxx, g_y_0_xzzzz_xxxxxy, g_y_0_xzzzz_xxxxxz, g_y_0_xzzzz_xxxxyy, g_y_0_xzzzz_xxxxyz, g_y_0_xzzzz_xxxxzz, g_y_0_xzzzz_xxxyyy, g_y_0_xzzzz_xxxyyz, g_y_0_xzzzz_xxxyzz, g_y_0_xzzzz_xxxzzz, g_y_0_xzzzz_xxyyyy, g_y_0_xzzzz_xxyyyz, g_y_0_xzzzz_xxyyzz, g_y_0_xzzzz_xxyzzz, g_y_0_xzzzz_xxzzzz, g_y_0_xzzzz_xyyyyy, g_y_0_xzzzz_xyyyyz, g_y_0_xzzzz_xyyyzz, g_y_0_xzzzz_xyyzzz, g_y_0_xzzzz_xyzzzz, g_y_0_xzzzz_xzzzzz, g_y_0_xzzzz_yyyyyy, g_y_0_xzzzz_yyyyyz, g_y_0_xzzzz_yyyyzz, g_y_0_xzzzz_yyyzzz, g_y_0_xzzzz_yyzzzz, g_y_0_xzzzz_yzzzzz, g_y_0_xzzzz_zzzzzz, g_y_0_zzzz_xxxxxx, g_y_0_zzzz_xxxxxxx, g_y_0_zzzz_xxxxxxy, g_y_0_zzzz_xxxxxxz, g_y_0_zzzz_xxxxxy, g_y_0_zzzz_xxxxxyy, g_y_0_zzzz_xxxxxyz, g_y_0_zzzz_xxxxxz, g_y_0_zzzz_xxxxxzz, g_y_0_zzzz_xxxxyy, g_y_0_zzzz_xxxxyyy, g_y_0_zzzz_xxxxyyz, g_y_0_zzzz_xxxxyz, g_y_0_zzzz_xxxxyzz, g_y_0_zzzz_xxxxzz, g_y_0_zzzz_xxxxzzz, g_y_0_zzzz_xxxyyy, g_y_0_zzzz_xxxyyyy, g_y_0_zzzz_xxxyyyz, g_y_0_zzzz_xxxyyz, g_y_0_zzzz_xxxyyzz, g_y_0_zzzz_xxxyzz, g_y_0_zzzz_xxxyzzz, g_y_0_zzzz_xxxzzz, g_y_0_zzzz_xxxzzzz, g_y_0_zzzz_xxyyyy, g_y_0_zzzz_xxyyyyy, g_y_0_zzzz_xxyyyyz, g_y_0_zzzz_xxyyyz, g_y_0_zzzz_xxyyyzz, g_y_0_zzzz_xxyyzz, g_y_0_zzzz_xxyyzzz, g_y_0_zzzz_xxyzzz, g_y_0_zzzz_xxyzzzz, g_y_0_zzzz_xxzzzz, g_y_0_zzzz_xxzzzzz, g_y_0_zzzz_xyyyyy, g_y_0_zzzz_xyyyyyy, g_y_0_zzzz_xyyyyyz, g_y_0_zzzz_xyyyyz, g_y_0_zzzz_xyyyyzz, g_y_0_zzzz_xyyyzz, g_y_0_zzzz_xyyyzzz, g_y_0_zzzz_xyyzzz, g_y_0_zzzz_xyyzzzz, g_y_0_zzzz_xyzzzz, g_y_0_zzzz_xyzzzzz, g_y_0_zzzz_xzzzzz, g_y_0_zzzz_xzzzzzz, g_y_0_zzzz_yyyyyy, g_y_0_zzzz_yyyyyz, g_y_0_zzzz_yyyyzz, g_y_0_zzzz_yyyzzz, g_y_0_zzzz_yyzzzz, g_y_0_zzzz_yzzzzz, g_y_0_zzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzz_xxxxxx[k] = -g_y_0_zzzz_xxxxxx[k] * cd_x[k] + g_y_0_zzzz_xxxxxxx[k];

                g_y_0_xzzzz_xxxxxy[k] = -g_y_0_zzzz_xxxxxy[k] * cd_x[k] + g_y_0_zzzz_xxxxxxy[k];

                g_y_0_xzzzz_xxxxxz[k] = -g_y_0_zzzz_xxxxxz[k] * cd_x[k] + g_y_0_zzzz_xxxxxxz[k];

                g_y_0_xzzzz_xxxxyy[k] = -g_y_0_zzzz_xxxxyy[k] * cd_x[k] + g_y_0_zzzz_xxxxxyy[k];

                g_y_0_xzzzz_xxxxyz[k] = -g_y_0_zzzz_xxxxyz[k] * cd_x[k] + g_y_0_zzzz_xxxxxyz[k];

                g_y_0_xzzzz_xxxxzz[k] = -g_y_0_zzzz_xxxxzz[k] * cd_x[k] + g_y_0_zzzz_xxxxxzz[k];

                g_y_0_xzzzz_xxxyyy[k] = -g_y_0_zzzz_xxxyyy[k] * cd_x[k] + g_y_0_zzzz_xxxxyyy[k];

                g_y_0_xzzzz_xxxyyz[k] = -g_y_0_zzzz_xxxyyz[k] * cd_x[k] + g_y_0_zzzz_xxxxyyz[k];

                g_y_0_xzzzz_xxxyzz[k] = -g_y_0_zzzz_xxxyzz[k] * cd_x[k] + g_y_0_zzzz_xxxxyzz[k];

                g_y_0_xzzzz_xxxzzz[k] = -g_y_0_zzzz_xxxzzz[k] * cd_x[k] + g_y_0_zzzz_xxxxzzz[k];

                g_y_0_xzzzz_xxyyyy[k] = -g_y_0_zzzz_xxyyyy[k] * cd_x[k] + g_y_0_zzzz_xxxyyyy[k];

                g_y_0_xzzzz_xxyyyz[k] = -g_y_0_zzzz_xxyyyz[k] * cd_x[k] + g_y_0_zzzz_xxxyyyz[k];

                g_y_0_xzzzz_xxyyzz[k] = -g_y_0_zzzz_xxyyzz[k] * cd_x[k] + g_y_0_zzzz_xxxyyzz[k];

                g_y_0_xzzzz_xxyzzz[k] = -g_y_0_zzzz_xxyzzz[k] * cd_x[k] + g_y_0_zzzz_xxxyzzz[k];

                g_y_0_xzzzz_xxzzzz[k] = -g_y_0_zzzz_xxzzzz[k] * cd_x[k] + g_y_0_zzzz_xxxzzzz[k];

                g_y_0_xzzzz_xyyyyy[k] = -g_y_0_zzzz_xyyyyy[k] * cd_x[k] + g_y_0_zzzz_xxyyyyy[k];

                g_y_0_xzzzz_xyyyyz[k] = -g_y_0_zzzz_xyyyyz[k] * cd_x[k] + g_y_0_zzzz_xxyyyyz[k];

                g_y_0_xzzzz_xyyyzz[k] = -g_y_0_zzzz_xyyyzz[k] * cd_x[k] + g_y_0_zzzz_xxyyyzz[k];

                g_y_0_xzzzz_xyyzzz[k] = -g_y_0_zzzz_xyyzzz[k] * cd_x[k] + g_y_0_zzzz_xxyyzzz[k];

                g_y_0_xzzzz_xyzzzz[k] = -g_y_0_zzzz_xyzzzz[k] * cd_x[k] + g_y_0_zzzz_xxyzzzz[k];

                g_y_0_xzzzz_xzzzzz[k] = -g_y_0_zzzz_xzzzzz[k] * cd_x[k] + g_y_0_zzzz_xxzzzzz[k];

                g_y_0_xzzzz_yyyyyy[k] = -g_y_0_zzzz_yyyyyy[k] * cd_x[k] + g_y_0_zzzz_xyyyyyy[k];

                g_y_0_xzzzz_yyyyyz[k] = -g_y_0_zzzz_yyyyyz[k] * cd_x[k] + g_y_0_zzzz_xyyyyyz[k];

                g_y_0_xzzzz_yyyyzz[k] = -g_y_0_zzzz_yyyyzz[k] * cd_x[k] + g_y_0_zzzz_xyyyyzz[k];

                g_y_0_xzzzz_yyyzzz[k] = -g_y_0_zzzz_yyyzzz[k] * cd_x[k] + g_y_0_zzzz_xyyyzzz[k];

                g_y_0_xzzzz_yyzzzz[k] = -g_y_0_zzzz_yyzzzz[k] * cd_x[k] + g_y_0_zzzz_xyyzzzz[k];

                g_y_0_xzzzz_yzzzzz[k] = -g_y_0_zzzz_yzzzzz[k] * cd_x[k] + g_y_0_zzzz_xyzzzzz[k];

                g_y_0_xzzzz_zzzzzz[k] = -g_y_0_zzzz_zzzzzz[k] * cd_x[k] + g_y_0_zzzz_xzzzzzz[k];
            }

            /// Set up 420-448 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 420);

            auto g_y_0_yyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 421);

            auto g_y_0_yyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 422);

            auto g_y_0_yyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 423);

            auto g_y_0_yyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 424);

            auto g_y_0_yyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 425);

            auto g_y_0_yyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 426);

            auto g_y_0_yyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 427);

            auto g_y_0_yyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 428);

            auto g_y_0_yyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 429);

            auto g_y_0_yyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 430);

            auto g_y_0_yyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 431);

            auto g_y_0_yyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 432);

            auto g_y_0_yyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 433);

            auto g_y_0_yyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 434);

            auto g_y_0_yyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 435);

            auto g_y_0_yyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 436);

            auto g_y_0_yyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 437);

            auto g_y_0_yyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 438);

            auto g_y_0_yyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 439);

            auto g_y_0_yyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 440);

            auto g_y_0_yyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 441);

            auto g_y_0_yyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 442);

            auto g_y_0_yyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 443);

            auto g_y_0_yyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 444);

            auto g_y_0_yyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 445);

            auto g_y_0_yyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 446);

            auto g_y_0_yyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 447);

            #pragma omp simd aligned(cd_y, g_y_0_yyyy_xxxxxx, g_y_0_yyyy_xxxxxxy, g_y_0_yyyy_xxxxxy, g_y_0_yyyy_xxxxxyy, g_y_0_yyyy_xxxxxyz, g_y_0_yyyy_xxxxxz, g_y_0_yyyy_xxxxyy, g_y_0_yyyy_xxxxyyy, g_y_0_yyyy_xxxxyyz, g_y_0_yyyy_xxxxyz, g_y_0_yyyy_xxxxyzz, g_y_0_yyyy_xxxxzz, g_y_0_yyyy_xxxyyy, g_y_0_yyyy_xxxyyyy, g_y_0_yyyy_xxxyyyz, g_y_0_yyyy_xxxyyz, g_y_0_yyyy_xxxyyzz, g_y_0_yyyy_xxxyzz, g_y_0_yyyy_xxxyzzz, g_y_0_yyyy_xxxzzz, g_y_0_yyyy_xxyyyy, g_y_0_yyyy_xxyyyyy, g_y_0_yyyy_xxyyyyz, g_y_0_yyyy_xxyyyz, g_y_0_yyyy_xxyyyzz, g_y_0_yyyy_xxyyzz, g_y_0_yyyy_xxyyzzz, g_y_0_yyyy_xxyzzz, g_y_0_yyyy_xxyzzzz, g_y_0_yyyy_xxzzzz, g_y_0_yyyy_xyyyyy, g_y_0_yyyy_xyyyyyy, g_y_0_yyyy_xyyyyyz, g_y_0_yyyy_xyyyyz, g_y_0_yyyy_xyyyyzz, g_y_0_yyyy_xyyyzz, g_y_0_yyyy_xyyyzzz, g_y_0_yyyy_xyyzzz, g_y_0_yyyy_xyyzzzz, g_y_0_yyyy_xyzzzz, g_y_0_yyyy_xyzzzzz, g_y_0_yyyy_xzzzzz, g_y_0_yyyy_yyyyyy, g_y_0_yyyy_yyyyyyy, g_y_0_yyyy_yyyyyyz, g_y_0_yyyy_yyyyyz, g_y_0_yyyy_yyyyyzz, g_y_0_yyyy_yyyyzz, g_y_0_yyyy_yyyyzzz, g_y_0_yyyy_yyyzzz, g_y_0_yyyy_yyyzzzz, g_y_0_yyyy_yyzzzz, g_y_0_yyyy_yyzzzzz, g_y_0_yyyy_yzzzzz, g_y_0_yyyy_yzzzzzz, g_y_0_yyyy_zzzzzz, g_y_0_yyyyy_xxxxxx, g_y_0_yyyyy_xxxxxy, g_y_0_yyyyy_xxxxxz, g_y_0_yyyyy_xxxxyy, g_y_0_yyyyy_xxxxyz, g_y_0_yyyyy_xxxxzz, g_y_0_yyyyy_xxxyyy, g_y_0_yyyyy_xxxyyz, g_y_0_yyyyy_xxxyzz, g_y_0_yyyyy_xxxzzz, g_y_0_yyyyy_xxyyyy, g_y_0_yyyyy_xxyyyz, g_y_0_yyyyy_xxyyzz, g_y_0_yyyyy_xxyzzz, g_y_0_yyyyy_xxzzzz, g_y_0_yyyyy_xyyyyy, g_y_0_yyyyy_xyyyyz, g_y_0_yyyyy_xyyyzz, g_y_0_yyyyy_xyyzzz, g_y_0_yyyyy_xyzzzz, g_y_0_yyyyy_xzzzzz, g_y_0_yyyyy_yyyyyy, g_y_0_yyyyy_yyyyyz, g_y_0_yyyyy_yyyyzz, g_y_0_yyyyy_yyyzzz, g_y_0_yyyyy_yyzzzz, g_y_0_yyyyy_yzzzzz, g_y_0_yyyyy_zzzzzz, g_yyyy_xxxxxx, g_yyyy_xxxxxy, g_yyyy_xxxxxz, g_yyyy_xxxxyy, g_yyyy_xxxxyz, g_yyyy_xxxxzz, g_yyyy_xxxyyy, g_yyyy_xxxyyz, g_yyyy_xxxyzz, g_yyyy_xxxzzz, g_yyyy_xxyyyy, g_yyyy_xxyyyz, g_yyyy_xxyyzz, g_yyyy_xxyzzz, g_yyyy_xxzzzz, g_yyyy_xyyyyy, g_yyyy_xyyyyz, g_yyyy_xyyyzz, g_yyyy_xyyzzz, g_yyyy_xyzzzz, g_yyyy_xzzzzz, g_yyyy_yyyyyy, g_yyyy_yyyyyz, g_yyyy_yyyyzz, g_yyyy_yyyzzz, g_yyyy_yyzzzz, g_yyyy_yzzzzz, g_yyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyy_xxxxxx[k] = -g_yyyy_xxxxxx[k] - g_y_0_yyyy_xxxxxx[k] * cd_y[k] + g_y_0_yyyy_xxxxxxy[k];

                g_y_0_yyyyy_xxxxxy[k] = -g_yyyy_xxxxxy[k] - g_y_0_yyyy_xxxxxy[k] * cd_y[k] + g_y_0_yyyy_xxxxxyy[k];

                g_y_0_yyyyy_xxxxxz[k] = -g_yyyy_xxxxxz[k] - g_y_0_yyyy_xxxxxz[k] * cd_y[k] + g_y_0_yyyy_xxxxxyz[k];

                g_y_0_yyyyy_xxxxyy[k] = -g_yyyy_xxxxyy[k] - g_y_0_yyyy_xxxxyy[k] * cd_y[k] + g_y_0_yyyy_xxxxyyy[k];

                g_y_0_yyyyy_xxxxyz[k] = -g_yyyy_xxxxyz[k] - g_y_0_yyyy_xxxxyz[k] * cd_y[k] + g_y_0_yyyy_xxxxyyz[k];

                g_y_0_yyyyy_xxxxzz[k] = -g_yyyy_xxxxzz[k] - g_y_0_yyyy_xxxxzz[k] * cd_y[k] + g_y_0_yyyy_xxxxyzz[k];

                g_y_0_yyyyy_xxxyyy[k] = -g_yyyy_xxxyyy[k] - g_y_0_yyyy_xxxyyy[k] * cd_y[k] + g_y_0_yyyy_xxxyyyy[k];

                g_y_0_yyyyy_xxxyyz[k] = -g_yyyy_xxxyyz[k] - g_y_0_yyyy_xxxyyz[k] * cd_y[k] + g_y_0_yyyy_xxxyyyz[k];

                g_y_0_yyyyy_xxxyzz[k] = -g_yyyy_xxxyzz[k] - g_y_0_yyyy_xxxyzz[k] * cd_y[k] + g_y_0_yyyy_xxxyyzz[k];

                g_y_0_yyyyy_xxxzzz[k] = -g_yyyy_xxxzzz[k] - g_y_0_yyyy_xxxzzz[k] * cd_y[k] + g_y_0_yyyy_xxxyzzz[k];

                g_y_0_yyyyy_xxyyyy[k] = -g_yyyy_xxyyyy[k] - g_y_0_yyyy_xxyyyy[k] * cd_y[k] + g_y_0_yyyy_xxyyyyy[k];

                g_y_0_yyyyy_xxyyyz[k] = -g_yyyy_xxyyyz[k] - g_y_0_yyyy_xxyyyz[k] * cd_y[k] + g_y_0_yyyy_xxyyyyz[k];

                g_y_0_yyyyy_xxyyzz[k] = -g_yyyy_xxyyzz[k] - g_y_0_yyyy_xxyyzz[k] * cd_y[k] + g_y_0_yyyy_xxyyyzz[k];

                g_y_0_yyyyy_xxyzzz[k] = -g_yyyy_xxyzzz[k] - g_y_0_yyyy_xxyzzz[k] * cd_y[k] + g_y_0_yyyy_xxyyzzz[k];

                g_y_0_yyyyy_xxzzzz[k] = -g_yyyy_xxzzzz[k] - g_y_0_yyyy_xxzzzz[k] * cd_y[k] + g_y_0_yyyy_xxyzzzz[k];

                g_y_0_yyyyy_xyyyyy[k] = -g_yyyy_xyyyyy[k] - g_y_0_yyyy_xyyyyy[k] * cd_y[k] + g_y_0_yyyy_xyyyyyy[k];

                g_y_0_yyyyy_xyyyyz[k] = -g_yyyy_xyyyyz[k] - g_y_0_yyyy_xyyyyz[k] * cd_y[k] + g_y_0_yyyy_xyyyyyz[k];

                g_y_0_yyyyy_xyyyzz[k] = -g_yyyy_xyyyzz[k] - g_y_0_yyyy_xyyyzz[k] * cd_y[k] + g_y_0_yyyy_xyyyyzz[k];

                g_y_0_yyyyy_xyyzzz[k] = -g_yyyy_xyyzzz[k] - g_y_0_yyyy_xyyzzz[k] * cd_y[k] + g_y_0_yyyy_xyyyzzz[k];

                g_y_0_yyyyy_xyzzzz[k] = -g_yyyy_xyzzzz[k] - g_y_0_yyyy_xyzzzz[k] * cd_y[k] + g_y_0_yyyy_xyyzzzz[k];

                g_y_0_yyyyy_xzzzzz[k] = -g_yyyy_xzzzzz[k] - g_y_0_yyyy_xzzzzz[k] * cd_y[k] + g_y_0_yyyy_xyzzzzz[k];

                g_y_0_yyyyy_yyyyyy[k] = -g_yyyy_yyyyyy[k] - g_y_0_yyyy_yyyyyy[k] * cd_y[k] + g_y_0_yyyy_yyyyyyy[k];

                g_y_0_yyyyy_yyyyyz[k] = -g_yyyy_yyyyyz[k] - g_y_0_yyyy_yyyyyz[k] * cd_y[k] + g_y_0_yyyy_yyyyyyz[k];

                g_y_0_yyyyy_yyyyzz[k] = -g_yyyy_yyyyzz[k] - g_y_0_yyyy_yyyyzz[k] * cd_y[k] + g_y_0_yyyy_yyyyyzz[k];

                g_y_0_yyyyy_yyyzzz[k] = -g_yyyy_yyyzzz[k] - g_y_0_yyyy_yyyzzz[k] * cd_y[k] + g_y_0_yyyy_yyyyzzz[k];

                g_y_0_yyyyy_yyzzzz[k] = -g_yyyy_yyzzzz[k] - g_y_0_yyyy_yyzzzz[k] * cd_y[k] + g_y_0_yyyy_yyyzzzz[k];

                g_y_0_yyyyy_yzzzzz[k] = -g_yyyy_yzzzzz[k] - g_y_0_yyyy_yzzzzz[k] * cd_y[k] + g_y_0_yyyy_yyzzzzz[k];

                g_y_0_yyyyy_zzzzzz[k] = -g_yyyy_zzzzzz[k] - g_y_0_yyyy_zzzzzz[k] * cd_y[k] + g_y_0_yyyy_yzzzzzz[k];
            }

            /// Set up 448-476 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 448);

            auto g_y_0_yyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 449);

            auto g_y_0_yyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 450);

            auto g_y_0_yyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 451);

            auto g_y_0_yyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 452);

            auto g_y_0_yyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 453);

            auto g_y_0_yyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 454);

            auto g_y_0_yyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 455);

            auto g_y_0_yyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 456);

            auto g_y_0_yyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 457);

            auto g_y_0_yyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 458);

            auto g_y_0_yyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 459);

            auto g_y_0_yyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 460);

            auto g_y_0_yyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 461);

            auto g_y_0_yyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 462);

            auto g_y_0_yyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 463);

            auto g_y_0_yyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 464);

            auto g_y_0_yyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 465);

            auto g_y_0_yyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 466);

            auto g_y_0_yyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 467);

            auto g_y_0_yyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 468);

            auto g_y_0_yyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 469);

            auto g_y_0_yyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 470);

            auto g_y_0_yyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 471);

            auto g_y_0_yyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 472);

            auto g_y_0_yyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 473);

            auto g_y_0_yyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 474);

            auto g_y_0_yyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 475);

            #pragma omp simd aligned(cd_z, g_y_0_yyyy_xxxxxx, g_y_0_yyyy_xxxxxxz, g_y_0_yyyy_xxxxxy, g_y_0_yyyy_xxxxxyz, g_y_0_yyyy_xxxxxz, g_y_0_yyyy_xxxxxzz, g_y_0_yyyy_xxxxyy, g_y_0_yyyy_xxxxyyz, g_y_0_yyyy_xxxxyz, g_y_0_yyyy_xxxxyzz, g_y_0_yyyy_xxxxzz, g_y_0_yyyy_xxxxzzz, g_y_0_yyyy_xxxyyy, g_y_0_yyyy_xxxyyyz, g_y_0_yyyy_xxxyyz, g_y_0_yyyy_xxxyyzz, g_y_0_yyyy_xxxyzz, g_y_0_yyyy_xxxyzzz, g_y_0_yyyy_xxxzzz, g_y_0_yyyy_xxxzzzz, g_y_0_yyyy_xxyyyy, g_y_0_yyyy_xxyyyyz, g_y_0_yyyy_xxyyyz, g_y_0_yyyy_xxyyyzz, g_y_0_yyyy_xxyyzz, g_y_0_yyyy_xxyyzzz, g_y_0_yyyy_xxyzzz, g_y_0_yyyy_xxyzzzz, g_y_0_yyyy_xxzzzz, g_y_0_yyyy_xxzzzzz, g_y_0_yyyy_xyyyyy, g_y_0_yyyy_xyyyyyz, g_y_0_yyyy_xyyyyz, g_y_0_yyyy_xyyyyzz, g_y_0_yyyy_xyyyzz, g_y_0_yyyy_xyyyzzz, g_y_0_yyyy_xyyzzz, g_y_0_yyyy_xyyzzzz, g_y_0_yyyy_xyzzzz, g_y_0_yyyy_xyzzzzz, g_y_0_yyyy_xzzzzz, g_y_0_yyyy_xzzzzzz, g_y_0_yyyy_yyyyyy, g_y_0_yyyy_yyyyyyz, g_y_0_yyyy_yyyyyz, g_y_0_yyyy_yyyyyzz, g_y_0_yyyy_yyyyzz, g_y_0_yyyy_yyyyzzz, g_y_0_yyyy_yyyzzz, g_y_0_yyyy_yyyzzzz, g_y_0_yyyy_yyzzzz, g_y_0_yyyy_yyzzzzz, g_y_0_yyyy_yzzzzz, g_y_0_yyyy_yzzzzzz, g_y_0_yyyy_zzzzzz, g_y_0_yyyy_zzzzzzz, g_y_0_yyyyz_xxxxxx, g_y_0_yyyyz_xxxxxy, g_y_0_yyyyz_xxxxxz, g_y_0_yyyyz_xxxxyy, g_y_0_yyyyz_xxxxyz, g_y_0_yyyyz_xxxxzz, g_y_0_yyyyz_xxxyyy, g_y_0_yyyyz_xxxyyz, g_y_0_yyyyz_xxxyzz, g_y_0_yyyyz_xxxzzz, g_y_0_yyyyz_xxyyyy, g_y_0_yyyyz_xxyyyz, g_y_0_yyyyz_xxyyzz, g_y_0_yyyyz_xxyzzz, g_y_0_yyyyz_xxzzzz, g_y_0_yyyyz_xyyyyy, g_y_0_yyyyz_xyyyyz, g_y_0_yyyyz_xyyyzz, g_y_0_yyyyz_xyyzzz, g_y_0_yyyyz_xyzzzz, g_y_0_yyyyz_xzzzzz, g_y_0_yyyyz_yyyyyy, g_y_0_yyyyz_yyyyyz, g_y_0_yyyyz_yyyyzz, g_y_0_yyyyz_yyyzzz, g_y_0_yyyyz_yyzzzz, g_y_0_yyyyz_yzzzzz, g_y_0_yyyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyz_xxxxxx[k] = -g_y_0_yyyy_xxxxxx[k] * cd_z[k] + g_y_0_yyyy_xxxxxxz[k];

                g_y_0_yyyyz_xxxxxy[k] = -g_y_0_yyyy_xxxxxy[k] * cd_z[k] + g_y_0_yyyy_xxxxxyz[k];

                g_y_0_yyyyz_xxxxxz[k] = -g_y_0_yyyy_xxxxxz[k] * cd_z[k] + g_y_0_yyyy_xxxxxzz[k];

                g_y_0_yyyyz_xxxxyy[k] = -g_y_0_yyyy_xxxxyy[k] * cd_z[k] + g_y_0_yyyy_xxxxyyz[k];

                g_y_0_yyyyz_xxxxyz[k] = -g_y_0_yyyy_xxxxyz[k] * cd_z[k] + g_y_0_yyyy_xxxxyzz[k];

                g_y_0_yyyyz_xxxxzz[k] = -g_y_0_yyyy_xxxxzz[k] * cd_z[k] + g_y_0_yyyy_xxxxzzz[k];

                g_y_0_yyyyz_xxxyyy[k] = -g_y_0_yyyy_xxxyyy[k] * cd_z[k] + g_y_0_yyyy_xxxyyyz[k];

                g_y_0_yyyyz_xxxyyz[k] = -g_y_0_yyyy_xxxyyz[k] * cd_z[k] + g_y_0_yyyy_xxxyyzz[k];

                g_y_0_yyyyz_xxxyzz[k] = -g_y_0_yyyy_xxxyzz[k] * cd_z[k] + g_y_0_yyyy_xxxyzzz[k];

                g_y_0_yyyyz_xxxzzz[k] = -g_y_0_yyyy_xxxzzz[k] * cd_z[k] + g_y_0_yyyy_xxxzzzz[k];

                g_y_0_yyyyz_xxyyyy[k] = -g_y_0_yyyy_xxyyyy[k] * cd_z[k] + g_y_0_yyyy_xxyyyyz[k];

                g_y_0_yyyyz_xxyyyz[k] = -g_y_0_yyyy_xxyyyz[k] * cd_z[k] + g_y_0_yyyy_xxyyyzz[k];

                g_y_0_yyyyz_xxyyzz[k] = -g_y_0_yyyy_xxyyzz[k] * cd_z[k] + g_y_0_yyyy_xxyyzzz[k];

                g_y_0_yyyyz_xxyzzz[k] = -g_y_0_yyyy_xxyzzz[k] * cd_z[k] + g_y_0_yyyy_xxyzzzz[k];

                g_y_0_yyyyz_xxzzzz[k] = -g_y_0_yyyy_xxzzzz[k] * cd_z[k] + g_y_0_yyyy_xxzzzzz[k];

                g_y_0_yyyyz_xyyyyy[k] = -g_y_0_yyyy_xyyyyy[k] * cd_z[k] + g_y_0_yyyy_xyyyyyz[k];

                g_y_0_yyyyz_xyyyyz[k] = -g_y_0_yyyy_xyyyyz[k] * cd_z[k] + g_y_0_yyyy_xyyyyzz[k];

                g_y_0_yyyyz_xyyyzz[k] = -g_y_0_yyyy_xyyyzz[k] * cd_z[k] + g_y_0_yyyy_xyyyzzz[k];

                g_y_0_yyyyz_xyyzzz[k] = -g_y_0_yyyy_xyyzzz[k] * cd_z[k] + g_y_0_yyyy_xyyzzzz[k];

                g_y_0_yyyyz_xyzzzz[k] = -g_y_0_yyyy_xyzzzz[k] * cd_z[k] + g_y_0_yyyy_xyzzzzz[k];

                g_y_0_yyyyz_xzzzzz[k] = -g_y_0_yyyy_xzzzzz[k] * cd_z[k] + g_y_0_yyyy_xzzzzzz[k];

                g_y_0_yyyyz_yyyyyy[k] = -g_y_0_yyyy_yyyyyy[k] * cd_z[k] + g_y_0_yyyy_yyyyyyz[k];

                g_y_0_yyyyz_yyyyyz[k] = -g_y_0_yyyy_yyyyyz[k] * cd_z[k] + g_y_0_yyyy_yyyyyzz[k];

                g_y_0_yyyyz_yyyyzz[k] = -g_y_0_yyyy_yyyyzz[k] * cd_z[k] + g_y_0_yyyy_yyyyzzz[k];

                g_y_0_yyyyz_yyyzzz[k] = -g_y_0_yyyy_yyyzzz[k] * cd_z[k] + g_y_0_yyyy_yyyzzzz[k];

                g_y_0_yyyyz_yyzzzz[k] = -g_y_0_yyyy_yyzzzz[k] * cd_z[k] + g_y_0_yyyy_yyzzzzz[k];

                g_y_0_yyyyz_yzzzzz[k] = -g_y_0_yyyy_yzzzzz[k] * cd_z[k] + g_y_0_yyyy_yzzzzzz[k];

                g_y_0_yyyyz_zzzzzz[k] = -g_y_0_yyyy_zzzzzz[k] * cd_z[k] + g_y_0_yyyy_zzzzzzz[k];
            }

            /// Set up 476-504 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 476);

            auto g_y_0_yyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 477);

            auto g_y_0_yyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 478);

            auto g_y_0_yyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 479);

            auto g_y_0_yyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 480);

            auto g_y_0_yyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 481);

            auto g_y_0_yyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 482);

            auto g_y_0_yyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 483);

            auto g_y_0_yyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 484);

            auto g_y_0_yyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 485);

            auto g_y_0_yyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 486);

            auto g_y_0_yyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 487);

            auto g_y_0_yyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 488);

            auto g_y_0_yyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 489);

            auto g_y_0_yyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 490);

            auto g_y_0_yyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 491);

            auto g_y_0_yyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 492);

            auto g_y_0_yyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 493);

            auto g_y_0_yyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 494);

            auto g_y_0_yyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 495);

            auto g_y_0_yyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 496);

            auto g_y_0_yyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 497);

            auto g_y_0_yyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 498);

            auto g_y_0_yyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 499);

            auto g_y_0_yyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 500);

            auto g_y_0_yyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 501);

            auto g_y_0_yyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 502);

            auto g_y_0_yyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 503);

            #pragma omp simd aligned(cd_z, g_y_0_yyyz_xxxxxx, g_y_0_yyyz_xxxxxxz, g_y_0_yyyz_xxxxxy, g_y_0_yyyz_xxxxxyz, g_y_0_yyyz_xxxxxz, g_y_0_yyyz_xxxxxzz, g_y_0_yyyz_xxxxyy, g_y_0_yyyz_xxxxyyz, g_y_0_yyyz_xxxxyz, g_y_0_yyyz_xxxxyzz, g_y_0_yyyz_xxxxzz, g_y_0_yyyz_xxxxzzz, g_y_0_yyyz_xxxyyy, g_y_0_yyyz_xxxyyyz, g_y_0_yyyz_xxxyyz, g_y_0_yyyz_xxxyyzz, g_y_0_yyyz_xxxyzz, g_y_0_yyyz_xxxyzzz, g_y_0_yyyz_xxxzzz, g_y_0_yyyz_xxxzzzz, g_y_0_yyyz_xxyyyy, g_y_0_yyyz_xxyyyyz, g_y_0_yyyz_xxyyyz, g_y_0_yyyz_xxyyyzz, g_y_0_yyyz_xxyyzz, g_y_0_yyyz_xxyyzzz, g_y_0_yyyz_xxyzzz, g_y_0_yyyz_xxyzzzz, g_y_0_yyyz_xxzzzz, g_y_0_yyyz_xxzzzzz, g_y_0_yyyz_xyyyyy, g_y_0_yyyz_xyyyyyz, g_y_0_yyyz_xyyyyz, g_y_0_yyyz_xyyyyzz, g_y_0_yyyz_xyyyzz, g_y_0_yyyz_xyyyzzz, g_y_0_yyyz_xyyzzz, g_y_0_yyyz_xyyzzzz, g_y_0_yyyz_xyzzzz, g_y_0_yyyz_xyzzzzz, g_y_0_yyyz_xzzzzz, g_y_0_yyyz_xzzzzzz, g_y_0_yyyz_yyyyyy, g_y_0_yyyz_yyyyyyz, g_y_0_yyyz_yyyyyz, g_y_0_yyyz_yyyyyzz, g_y_0_yyyz_yyyyzz, g_y_0_yyyz_yyyyzzz, g_y_0_yyyz_yyyzzz, g_y_0_yyyz_yyyzzzz, g_y_0_yyyz_yyzzzz, g_y_0_yyyz_yyzzzzz, g_y_0_yyyz_yzzzzz, g_y_0_yyyz_yzzzzzz, g_y_0_yyyz_zzzzzz, g_y_0_yyyz_zzzzzzz, g_y_0_yyyzz_xxxxxx, g_y_0_yyyzz_xxxxxy, g_y_0_yyyzz_xxxxxz, g_y_0_yyyzz_xxxxyy, g_y_0_yyyzz_xxxxyz, g_y_0_yyyzz_xxxxzz, g_y_0_yyyzz_xxxyyy, g_y_0_yyyzz_xxxyyz, g_y_0_yyyzz_xxxyzz, g_y_0_yyyzz_xxxzzz, g_y_0_yyyzz_xxyyyy, g_y_0_yyyzz_xxyyyz, g_y_0_yyyzz_xxyyzz, g_y_0_yyyzz_xxyzzz, g_y_0_yyyzz_xxzzzz, g_y_0_yyyzz_xyyyyy, g_y_0_yyyzz_xyyyyz, g_y_0_yyyzz_xyyyzz, g_y_0_yyyzz_xyyzzz, g_y_0_yyyzz_xyzzzz, g_y_0_yyyzz_xzzzzz, g_y_0_yyyzz_yyyyyy, g_y_0_yyyzz_yyyyyz, g_y_0_yyyzz_yyyyzz, g_y_0_yyyzz_yyyzzz, g_y_0_yyyzz_yyzzzz, g_y_0_yyyzz_yzzzzz, g_y_0_yyyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzz_xxxxxx[k] = -g_y_0_yyyz_xxxxxx[k] * cd_z[k] + g_y_0_yyyz_xxxxxxz[k];

                g_y_0_yyyzz_xxxxxy[k] = -g_y_0_yyyz_xxxxxy[k] * cd_z[k] + g_y_0_yyyz_xxxxxyz[k];

                g_y_0_yyyzz_xxxxxz[k] = -g_y_0_yyyz_xxxxxz[k] * cd_z[k] + g_y_0_yyyz_xxxxxzz[k];

                g_y_0_yyyzz_xxxxyy[k] = -g_y_0_yyyz_xxxxyy[k] * cd_z[k] + g_y_0_yyyz_xxxxyyz[k];

                g_y_0_yyyzz_xxxxyz[k] = -g_y_0_yyyz_xxxxyz[k] * cd_z[k] + g_y_0_yyyz_xxxxyzz[k];

                g_y_0_yyyzz_xxxxzz[k] = -g_y_0_yyyz_xxxxzz[k] * cd_z[k] + g_y_0_yyyz_xxxxzzz[k];

                g_y_0_yyyzz_xxxyyy[k] = -g_y_0_yyyz_xxxyyy[k] * cd_z[k] + g_y_0_yyyz_xxxyyyz[k];

                g_y_0_yyyzz_xxxyyz[k] = -g_y_0_yyyz_xxxyyz[k] * cd_z[k] + g_y_0_yyyz_xxxyyzz[k];

                g_y_0_yyyzz_xxxyzz[k] = -g_y_0_yyyz_xxxyzz[k] * cd_z[k] + g_y_0_yyyz_xxxyzzz[k];

                g_y_0_yyyzz_xxxzzz[k] = -g_y_0_yyyz_xxxzzz[k] * cd_z[k] + g_y_0_yyyz_xxxzzzz[k];

                g_y_0_yyyzz_xxyyyy[k] = -g_y_0_yyyz_xxyyyy[k] * cd_z[k] + g_y_0_yyyz_xxyyyyz[k];

                g_y_0_yyyzz_xxyyyz[k] = -g_y_0_yyyz_xxyyyz[k] * cd_z[k] + g_y_0_yyyz_xxyyyzz[k];

                g_y_0_yyyzz_xxyyzz[k] = -g_y_0_yyyz_xxyyzz[k] * cd_z[k] + g_y_0_yyyz_xxyyzzz[k];

                g_y_0_yyyzz_xxyzzz[k] = -g_y_0_yyyz_xxyzzz[k] * cd_z[k] + g_y_0_yyyz_xxyzzzz[k];

                g_y_0_yyyzz_xxzzzz[k] = -g_y_0_yyyz_xxzzzz[k] * cd_z[k] + g_y_0_yyyz_xxzzzzz[k];

                g_y_0_yyyzz_xyyyyy[k] = -g_y_0_yyyz_xyyyyy[k] * cd_z[k] + g_y_0_yyyz_xyyyyyz[k];

                g_y_0_yyyzz_xyyyyz[k] = -g_y_0_yyyz_xyyyyz[k] * cd_z[k] + g_y_0_yyyz_xyyyyzz[k];

                g_y_0_yyyzz_xyyyzz[k] = -g_y_0_yyyz_xyyyzz[k] * cd_z[k] + g_y_0_yyyz_xyyyzzz[k];

                g_y_0_yyyzz_xyyzzz[k] = -g_y_0_yyyz_xyyzzz[k] * cd_z[k] + g_y_0_yyyz_xyyzzzz[k];

                g_y_0_yyyzz_xyzzzz[k] = -g_y_0_yyyz_xyzzzz[k] * cd_z[k] + g_y_0_yyyz_xyzzzzz[k];

                g_y_0_yyyzz_xzzzzz[k] = -g_y_0_yyyz_xzzzzz[k] * cd_z[k] + g_y_0_yyyz_xzzzzzz[k];

                g_y_0_yyyzz_yyyyyy[k] = -g_y_0_yyyz_yyyyyy[k] * cd_z[k] + g_y_0_yyyz_yyyyyyz[k];

                g_y_0_yyyzz_yyyyyz[k] = -g_y_0_yyyz_yyyyyz[k] * cd_z[k] + g_y_0_yyyz_yyyyyzz[k];

                g_y_0_yyyzz_yyyyzz[k] = -g_y_0_yyyz_yyyyzz[k] * cd_z[k] + g_y_0_yyyz_yyyyzzz[k];

                g_y_0_yyyzz_yyyzzz[k] = -g_y_0_yyyz_yyyzzz[k] * cd_z[k] + g_y_0_yyyz_yyyzzzz[k];

                g_y_0_yyyzz_yyzzzz[k] = -g_y_0_yyyz_yyzzzz[k] * cd_z[k] + g_y_0_yyyz_yyzzzzz[k];

                g_y_0_yyyzz_yzzzzz[k] = -g_y_0_yyyz_yzzzzz[k] * cd_z[k] + g_y_0_yyyz_yzzzzzz[k];

                g_y_0_yyyzz_zzzzzz[k] = -g_y_0_yyyz_zzzzzz[k] * cd_z[k] + g_y_0_yyyz_zzzzzzz[k];
            }

            /// Set up 504-532 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 504);

            auto g_y_0_yyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 505);

            auto g_y_0_yyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 506);

            auto g_y_0_yyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 507);

            auto g_y_0_yyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 508);

            auto g_y_0_yyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 509);

            auto g_y_0_yyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 510);

            auto g_y_0_yyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 511);

            auto g_y_0_yyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 512);

            auto g_y_0_yyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 513);

            auto g_y_0_yyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 514);

            auto g_y_0_yyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 515);

            auto g_y_0_yyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 516);

            auto g_y_0_yyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 517);

            auto g_y_0_yyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 518);

            auto g_y_0_yyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 519);

            auto g_y_0_yyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 520);

            auto g_y_0_yyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 521);

            auto g_y_0_yyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 522);

            auto g_y_0_yyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 523);

            auto g_y_0_yyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 524);

            auto g_y_0_yyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 525);

            auto g_y_0_yyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 526);

            auto g_y_0_yyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 527);

            auto g_y_0_yyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 528);

            auto g_y_0_yyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 529);

            auto g_y_0_yyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 530);

            auto g_y_0_yyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 531);

            #pragma omp simd aligned(cd_z, g_y_0_yyzz_xxxxxx, g_y_0_yyzz_xxxxxxz, g_y_0_yyzz_xxxxxy, g_y_0_yyzz_xxxxxyz, g_y_0_yyzz_xxxxxz, g_y_0_yyzz_xxxxxzz, g_y_0_yyzz_xxxxyy, g_y_0_yyzz_xxxxyyz, g_y_0_yyzz_xxxxyz, g_y_0_yyzz_xxxxyzz, g_y_0_yyzz_xxxxzz, g_y_0_yyzz_xxxxzzz, g_y_0_yyzz_xxxyyy, g_y_0_yyzz_xxxyyyz, g_y_0_yyzz_xxxyyz, g_y_0_yyzz_xxxyyzz, g_y_0_yyzz_xxxyzz, g_y_0_yyzz_xxxyzzz, g_y_0_yyzz_xxxzzz, g_y_0_yyzz_xxxzzzz, g_y_0_yyzz_xxyyyy, g_y_0_yyzz_xxyyyyz, g_y_0_yyzz_xxyyyz, g_y_0_yyzz_xxyyyzz, g_y_0_yyzz_xxyyzz, g_y_0_yyzz_xxyyzzz, g_y_0_yyzz_xxyzzz, g_y_0_yyzz_xxyzzzz, g_y_0_yyzz_xxzzzz, g_y_0_yyzz_xxzzzzz, g_y_0_yyzz_xyyyyy, g_y_0_yyzz_xyyyyyz, g_y_0_yyzz_xyyyyz, g_y_0_yyzz_xyyyyzz, g_y_0_yyzz_xyyyzz, g_y_0_yyzz_xyyyzzz, g_y_0_yyzz_xyyzzz, g_y_0_yyzz_xyyzzzz, g_y_0_yyzz_xyzzzz, g_y_0_yyzz_xyzzzzz, g_y_0_yyzz_xzzzzz, g_y_0_yyzz_xzzzzzz, g_y_0_yyzz_yyyyyy, g_y_0_yyzz_yyyyyyz, g_y_0_yyzz_yyyyyz, g_y_0_yyzz_yyyyyzz, g_y_0_yyzz_yyyyzz, g_y_0_yyzz_yyyyzzz, g_y_0_yyzz_yyyzzz, g_y_0_yyzz_yyyzzzz, g_y_0_yyzz_yyzzzz, g_y_0_yyzz_yyzzzzz, g_y_0_yyzz_yzzzzz, g_y_0_yyzz_yzzzzzz, g_y_0_yyzz_zzzzzz, g_y_0_yyzz_zzzzzzz, g_y_0_yyzzz_xxxxxx, g_y_0_yyzzz_xxxxxy, g_y_0_yyzzz_xxxxxz, g_y_0_yyzzz_xxxxyy, g_y_0_yyzzz_xxxxyz, g_y_0_yyzzz_xxxxzz, g_y_0_yyzzz_xxxyyy, g_y_0_yyzzz_xxxyyz, g_y_0_yyzzz_xxxyzz, g_y_0_yyzzz_xxxzzz, g_y_0_yyzzz_xxyyyy, g_y_0_yyzzz_xxyyyz, g_y_0_yyzzz_xxyyzz, g_y_0_yyzzz_xxyzzz, g_y_0_yyzzz_xxzzzz, g_y_0_yyzzz_xyyyyy, g_y_0_yyzzz_xyyyyz, g_y_0_yyzzz_xyyyzz, g_y_0_yyzzz_xyyzzz, g_y_0_yyzzz_xyzzzz, g_y_0_yyzzz_xzzzzz, g_y_0_yyzzz_yyyyyy, g_y_0_yyzzz_yyyyyz, g_y_0_yyzzz_yyyyzz, g_y_0_yyzzz_yyyzzz, g_y_0_yyzzz_yyzzzz, g_y_0_yyzzz_yzzzzz, g_y_0_yyzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzz_xxxxxx[k] = -g_y_0_yyzz_xxxxxx[k] * cd_z[k] + g_y_0_yyzz_xxxxxxz[k];

                g_y_0_yyzzz_xxxxxy[k] = -g_y_0_yyzz_xxxxxy[k] * cd_z[k] + g_y_0_yyzz_xxxxxyz[k];

                g_y_0_yyzzz_xxxxxz[k] = -g_y_0_yyzz_xxxxxz[k] * cd_z[k] + g_y_0_yyzz_xxxxxzz[k];

                g_y_0_yyzzz_xxxxyy[k] = -g_y_0_yyzz_xxxxyy[k] * cd_z[k] + g_y_0_yyzz_xxxxyyz[k];

                g_y_0_yyzzz_xxxxyz[k] = -g_y_0_yyzz_xxxxyz[k] * cd_z[k] + g_y_0_yyzz_xxxxyzz[k];

                g_y_0_yyzzz_xxxxzz[k] = -g_y_0_yyzz_xxxxzz[k] * cd_z[k] + g_y_0_yyzz_xxxxzzz[k];

                g_y_0_yyzzz_xxxyyy[k] = -g_y_0_yyzz_xxxyyy[k] * cd_z[k] + g_y_0_yyzz_xxxyyyz[k];

                g_y_0_yyzzz_xxxyyz[k] = -g_y_0_yyzz_xxxyyz[k] * cd_z[k] + g_y_0_yyzz_xxxyyzz[k];

                g_y_0_yyzzz_xxxyzz[k] = -g_y_0_yyzz_xxxyzz[k] * cd_z[k] + g_y_0_yyzz_xxxyzzz[k];

                g_y_0_yyzzz_xxxzzz[k] = -g_y_0_yyzz_xxxzzz[k] * cd_z[k] + g_y_0_yyzz_xxxzzzz[k];

                g_y_0_yyzzz_xxyyyy[k] = -g_y_0_yyzz_xxyyyy[k] * cd_z[k] + g_y_0_yyzz_xxyyyyz[k];

                g_y_0_yyzzz_xxyyyz[k] = -g_y_0_yyzz_xxyyyz[k] * cd_z[k] + g_y_0_yyzz_xxyyyzz[k];

                g_y_0_yyzzz_xxyyzz[k] = -g_y_0_yyzz_xxyyzz[k] * cd_z[k] + g_y_0_yyzz_xxyyzzz[k];

                g_y_0_yyzzz_xxyzzz[k] = -g_y_0_yyzz_xxyzzz[k] * cd_z[k] + g_y_0_yyzz_xxyzzzz[k];

                g_y_0_yyzzz_xxzzzz[k] = -g_y_0_yyzz_xxzzzz[k] * cd_z[k] + g_y_0_yyzz_xxzzzzz[k];

                g_y_0_yyzzz_xyyyyy[k] = -g_y_0_yyzz_xyyyyy[k] * cd_z[k] + g_y_0_yyzz_xyyyyyz[k];

                g_y_0_yyzzz_xyyyyz[k] = -g_y_0_yyzz_xyyyyz[k] * cd_z[k] + g_y_0_yyzz_xyyyyzz[k];

                g_y_0_yyzzz_xyyyzz[k] = -g_y_0_yyzz_xyyyzz[k] * cd_z[k] + g_y_0_yyzz_xyyyzzz[k];

                g_y_0_yyzzz_xyyzzz[k] = -g_y_0_yyzz_xyyzzz[k] * cd_z[k] + g_y_0_yyzz_xyyzzzz[k];

                g_y_0_yyzzz_xyzzzz[k] = -g_y_0_yyzz_xyzzzz[k] * cd_z[k] + g_y_0_yyzz_xyzzzzz[k];

                g_y_0_yyzzz_xzzzzz[k] = -g_y_0_yyzz_xzzzzz[k] * cd_z[k] + g_y_0_yyzz_xzzzzzz[k];

                g_y_0_yyzzz_yyyyyy[k] = -g_y_0_yyzz_yyyyyy[k] * cd_z[k] + g_y_0_yyzz_yyyyyyz[k];

                g_y_0_yyzzz_yyyyyz[k] = -g_y_0_yyzz_yyyyyz[k] * cd_z[k] + g_y_0_yyzz_yyyyyzz[k];

                g_y_0_yyzzz_yyyyzz[k] = -g_y_0_yyzz_yyyyzz[k] * cd_z[k] + g_y_0_yyzz_yyyyzzz[k];

                g_y_0_yyzzz_yyyzzz[k] = -g_y_0_yyzz_yyyzzz[k] * cd_z[k] + g_y_0_yyzz_yyyzzzz[k];

                g_y_0_yyzzz_yyzzzz[k] = -g_y_0_yyzz_yyzzzz[k] * cd_z[k] + g_y_0_yyzz_yyzzzzz[k];

                g_y_0_yyzzz_yzzzzz[k] = -g_y_0_yyzz_yzzzzz[k] * cd_z[k] + g_y_0_yyzz_yzzzzzz[k];

                g_y_0_yyzzz_zzzzzz[k] = -g_y_0_yyzz_zzzzzz[k] * cd_z[k] + g_y_0_yyzz_zzzzzzz[k];
            }

            /// Set up 532-560 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 532);

            auto g_y_0_yzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 533);

            auto g_y_0_yzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 534);

            auto g_y_0_yzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 535);

            auto g_y_0_yzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 536);

            auto g_y_0_yzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 537);

            auto g_y_0_yzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 538);

            auto g_y_0_yzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 539);

            auto g_y_0_yzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 540);

            auto g_y_0_yzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 541);

            auto g_y_0_yzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 542);

            auto g_y_0_yzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 543);

            auto g_y_0_yzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 544);

            auto g_y_0_yzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 545);

            auto g_y_0_yzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 546);

            auto g_y_0_yzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 547);

            auto g_y_0_yzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 548);

            auto g_y_0_yzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 549);

            auto g_y_0_yzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 550);

            auto g_y_0_yzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 551);

            auto g_y_0_yzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 552);

            auto g_y_0_yzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 553);

            auto g_y_0_yzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 554);

            auto g_y_0_yzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 555);

            auto g_y_0_yzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 556);

            auto g_y_0_yzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 557);

            auto g_y_0_yzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 558);

            auto g_y_0_yzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 559);

            #pragma omp simd aligned(cd_z, g_y_0_yzzz_xxxxxx, g_y_0_yzzz_xxxxxxz, g_y_0_yzzz_xxxxxy, g_y_0_yzzz_xxxxxyz, g_y_0_yzzz_xxxxxz, g_y_0_yzzz_xxxxxzz, g_y_0_yzzz_xxxxyy, g_y_0_yzzz_xxxxyyz, g_y_0_yzzz_xxxxyz, g_y_0_yzzz_xxxxyzz, g_y_0_yzzz_xxxxzz, g_y_0_yzzz_xxxxzzz, g_y_0_yzzz_xxxyyy, g_y_0_yzzz_xxxyyyz, g_y_0_yzzz_xxxyyz, g_y_0_yzzz_xxxyyzz, g_y_0_yzzz_xxxyzz, g_y_0_yzzz_xxxyzzz, g_y_0_yzzz_xxxzzz, g_y_0_yzzz_xxxzzzz, g_y_0_yzzz_xxyyyy, g_y_0_yzzz_xxyyyyz, g_y_0_yzzz_xxyyyz, g_y_0_yzzz_xxyyyzz, g_y_0_yzzz_xxyyzz, g_y_0_yzzz_xxyyzzz, g_y_0_yzzz_xxyzzz, g_y_0_yzzz_xxyzzzz, g_y_0_yzzz_xxzzzz, g_y_0_yzzz_xxzzzzz, g_y_0_yzzz_xyyyyy, g_y_0_yzzz_xyyyyyz, g_y_0_yzzz_xyyyyz, g_y_0_yzzz_xyyyyzz, g_y_0_yzzz_xyyyzz, g_y_0_yzzz_xyyyzzz, g_y_0_yzzz_xyyzzz, g_y_0_yzzz_xyyzzzz, g_y_0_yzzz_xyzzzz, g_y_0_yzzz_xyzzzzz, g_y_0_yzzz_xzzzzz, g_y_0_yzzz_xzzzzzz, g_y_0_yzzz_yyyyyy, g_y_0_yzzz_yyyyyyz, g_y_0_yzzz_yyyyyz, g_y_0_yzzz_yyyyyzz, g_y_0_yzzz_yyyyzz, g_y_0_yzzz_yyyyzzz, g_y_0_yzzz_yyyzzz, g_y_0_yzzz_yyyzzzz, g_y_0_yzzz_yyzzzz, g_y_0_yzzz_yyzzzzz, g_y_0_yzzz_yzzzzz, g_y_0_yzzz_yzzzzzz, g_y_0_yzzz_zzzzzz, g_y_0_yzzz_zzzzzzz, g_y_0_yzzzz_xxxxxx, g_y_0_yzzzz_xxxxxy, g_y_0_yzzzz_xxxxxz, g_y_0_yzzzz_xxxxyy, g_y_0_yzzzz_xxxxyz, g_y_0_yzzzz_xxxxzz, g_y_0_yzzzz_xxxyyy, g_y_0_yzzzz_xxxyyz, g_y_0_yzzzz_xxxyzz, g_y_0_yzzzz_xxxzzz, g_y_0_yzzzz_xxyyyy, g_y_0_yzzzz_xxyyyz, g_y_0_yzzzz_xxyyzz, g_y_0_yzzzz_xxyzzz, g_y_0_yzzzz_xxzzzz, g_y_0_yzzzz_xyyyyy, g_y_0_yzzzz_xyyyyz, g_y_0_yzzzz_xyyyzz, g_y_0_yzzzz_xyyzzz, g_y_0_yzzzz_xyzzzz, g_y_0_yzzzz_xzzzzz, g_y_0_yzzzz_yyyyyy, g_y_0_yzzzz_yyyyyz, g_y_0_yzzzz_yyyyzz, g_y_0_yzzzz_yyyzzz, g_y_0_yzzzz_yyzzzz, g_y_0_yzzzz_yzzzzz, g_y_0_yzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzz_xxxxxx[k] = -g_y_0_yzzz_xxxxxx[k] * cd_z[k] + g_y_0_yzzz_xxxxxxz[k];

                g_y_0_yzzzz_xxxxxy[k] = -g_y_0_yzzz_xxxxxy[k] * cd_z[k] + g_y_0_yzzz_xxxxxyz[k];

                g_y_0_yzzzz_xxxxxz[k] = -g_y_0_yzzz_xxxxxz[k] * cd_z[k] + g_y_0_yzzz_xxxxxzz[k];

                g_y_0_yzzzz_xxxxyy[k] = -g_y_0_yzzz_xxxxyy[k] * cd_z[k] + g_y_0_yzzz_xxxxyyz[k];

                g_y_0_yzzzz_xxxxyz[k] = -g_y_0_yzzz_xxxxyz[k] * cd_z[k] + g_y_0_yzzz_xxxxyzz[k];

                g_y_0_yzzzz_xxxxzz[k] = -g_y_0_yzzz_xxxxzz[k] * cd_z[k] + g_y_0_yzzz_xxxxzzz[k];

                g_y_0_yzzzz_xxxyyy[k] = -g_y_0_yzzz_xxxyyy[k] * cd_z[k] + g_y_0_yzzz_xxxyyyz[k];

                g_y_0_yzzzz_xxxyyz[k] = -g_y_0_yzzz_xxxyyz[k] * cd_z[k] + g_y_0_yzzz_xxxyyzz[k];

                g_y_0_yzzzz_xxxyzz[k] = -g_y_0_yzzz_xxxyzz[k] * cd_z[k] + g_y_0_yzzz_xxxyzzz[k];

                g_y_0_yzzzz_xxxzzz[k] = -g_y_0_yzzz_xxxzzz[k] * cd_z[k] + g_y_0_yzzz_xxxzzzz[k];

                g_y_0_yzzzz_xxyyyy[k] = -g_y_0_yzzz_xxyyyy[k] * cd_z[k] + g_y_0_yzzz_xxyyyyz[k];

                g_y_0_yzzzz_xxyyyz[k] = -g_y_0_yzzz_xxyyyz[k] * cd_z[k] + g_y_0_yzzz_xxyyyzz[k];

                g_y_0_yzzzz_xxyyzz[k] = -g_y_0_yzzz_xxyyzz[k] * cd_z[k] + g_y_0_yzzz_xxyyzzz[k];

                g_y_0_yzzzz_xxyzzz[k] = -g_y_0_yzzz_xxyzzz[k] * cd_z[k] + g_y_0_yzzz_xxyzzzz[k];

                g_y_0_yzzzz_xxzzzz[k] = -g_y_0_yzzz_xxzzzz[k] * cd_z[k] + g_y_0_yzzz_xxzzzzz[k];

                g_y_0_yzzzz_xyyyyy[k] = -g_y_0_yzzz_xyyyyy[k] * cd_z[k] + g_y_0_yzzz_xyyyyyz[k];

                g_y_0_yzzzz_xyyyyz[k] = -g_y_0_yzzz_xyyyyz[k] * cd_z[k] + g_y_0_yzzz_xyyyyzz[k];

                g_y_0_yzzzz_xyyyzz[k] = -g_y_0_yzzz_xyyyzz[k] * cd_z[k] + g_y_0_yzzz_xyyyzzz[k];

                g_y_0_yzzzz_xyyzzz[k] = -g_y_0_yzzz_xyyzzz[k] * cd_z[k] + g_y_0_yzzz_xyyzzzz[k];

                g_y_0_yzzzz_xyzzzz[k] = -g_y_0_yzzz_xyzzzz[k] * cd_z[k] + g_y_0_yzzz_xyzzzzz[k];

                g_y_0_yzzzz_xzzzzz[k] = -g_y_0_yzzz_xzzzzz[k] * cd_z[k] + g_y_0_yzzz_xzzzzzz[k];

                g_y_0_yzzzz_yyyyyy[k] = -g_y_0_yzzz_yyyyyy[k] * cd_z[k] + g_y_0_yzzz_yyyyyyz[k];

                g_y_0_yzzzz_yyyyyz[k] = -g_y_0_yzzz_yyyyyz[k] * cd_z[k] + g_y_0_yzzz_yyyyyzz[k];

                g_y_0_yzzzz_yyyyzz[k] = -g_y_0_yzzz_yyyyzz[k] * cd_z[k] + g_y_0_yzzz_yyyyzzz[k];

                g_y_0_yzzzz_yyyzzz[k] = -g_y_0_yzzz_yyyzzz[k] * cd_z[k] + g_y_0_yzzz_yyyzzzz[k];

                g_y_0_yzzzz_yyzzzz[k] = -g_y_0_yzzz_yyzzzz[k] * cd_z[k] + g_y_0_yzzz_yyzzzzz[k];

                g_y_0_yzzzz_yzzzzz[k] = -g_y_0_yzzz_yzzzzz[k] * cd_z[k] + g_y_0_yzzz_yzzzzzz[k];

                g_y_0_yzzzz_zzzzzz[k] = -g_y_0_yzzz_zzzzzz[k] * cd_z[k] + g_y_0_yzzz_zzzzzzz[k];
            }

            /// Set up 560-588 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 560);

            auto g_y_0_zzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 561);

            auto g_y_0_zzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 562);

            auto g_y_0_zzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 563);

            auto g_y_0_zzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 564);

            auto g_y_0_zzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 565);

            auto g_y_0_zzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 566);

            auto g_y_0_zzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 567);

            auto g_y_0_zzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 568);

            auto g_y_0_zzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 569);

            auto g_y_0_zzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 570);

            auto g_y_0_zzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 571);

            auto g_y_0_zzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 572);

            auto g_y_0_zzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 573);

            auto g_y_0_zzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 574);

            auto g_y_0_zzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 575);

            auto g_y_0_zzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 576);

            auto g_y_0_zzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 577);

            auto g_y_0_zzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 578);

            auto g_y_0_zzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 579);

            auto g_y_0_zzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 580);

            auto g_y_0_zzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 581);

            auto g_y_0_zzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 582);

            auto g_y_0_zzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 583);

            auto g_y_0_zzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 584);

            auto g_y_0_zzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 585);

            auto g_y_0_zzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 586);

            auto g_y_0_zzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 587);

            #pragma omp simd aligned(cd_z, g_y_0_zzzz_xxxxxx, g_y_0_zzzz_xxxxxxz, g_y_0_zzzz_xxxxxy, g_y_0_zzzz_xxxxxyz, g_y_0_zzzz_xxxxxz, g_y_0_zzzz_xxxxxzz, g_y_0_zzzz_xxxxyy, g_y_0_zzzz_xxxxyyz, g_y_0_zzzz_xxxxyz, g_y_0_zzzz_xxxxyzz, g_y_0_zzzz_xxxxzz, g_y_0_zzzz_xxxxzzz, g_y_0_zzzz_xxxyyy, g_y_0_zzzz_xxxyyyz, g_y_0_zzzz_xxxyyz, g_y_0_zzzz_xxxyyzz, g_y_0_zzzz_xxxyzz, g_y_0_zzzz_xxxyzzz, g_y_0_zzzz_xxxzzz, g_y_0_zzzz_xxxzzzz, g_y_0_zzzz_xxyyyy, g_y_0_zzzz_xxyyyyz, g_y_0_zzzz_xxyyyz, g_y_0_zzzz_xxyyyzz, g_y_0_zzzz_xxyyzz, g_y_0_zzzz_xxyyzzz, g_y_0_zzzz_xxyzzz, g_y_0_zzzz_xxyzzzz, g_y_0_zzzz_xxzzzz, g_y_0_zzzz_xxzzzzz, g_y_0_zzzz_xyyyyy, g_y_0_zzzz_xyyyyyz, g_y_0_zzzz_xyyyyz, g_y_0_zzzz_xyyyyzz, g_y_0_zzzz_xyyyzz, g_y_0_zzzz_xyyyzzz, g_y_0_zzzz_xyyzzz, g_y_0_zzzz_xyyzzzz, g_y_0_zzzz_xyzzzz, g_y_0_zzzz_xyzzzzz, g_y_0_zzzz_xzzzzz, g_y_0_zzzz_xzzzzzz, g_y_0_zzzz_yyyyyy, g_y_0_zzzz_yyyyyyz, g_y_0_zzzz_yyyyyz, g_y_0_zzzz_yyyyyzz, g_y_0_zzzz_yyyyzz, g_y_0_zzzz_yyyyzzz, g_y_0_zzzz_yyyzzz, g_y_0_zzzz_yyyzzzz, g_y_0_zzzz_yyzzzz, g_y_0_zzzz_yyzzzzz, g_y_0_zzzz_yzzzzz, g_y_0_zzzz_yzzzzzz, g_y_0_zzzz_zzzzzz, g_y_0_zzzz_zzzzzzz, g_y_0_zzzzz_xxxxxx, g_y_0_zzzzz_xxxxxy, g_y_0_zzzzz_xxxxxz, g_y_0_zzzzz_xxxxyy, g_y_0_zzzzz_xxxxyz, g_y_0_zzzzz_xxxxzz, g_y_0_zzzzz_xxxyyy, g_y_0_zzzzz_xxxyyz, g_y_0_zzzzz_xxxyzz, g_y_0_zzzzz_xxxzzz, g_y_0_zzzzz_xxyyyy, g_y_0_zzzzz_xxyyyz, g_y_0_zzzzz_xxyyzz, g_y_0_zzzzz_xxyzzz, g_y_0_zzzzz_xxzzzz, g_y_0_zzzzz_xyyyyy, g_y_0_zzzzz_xyyyyz, g_y_0_zzzzz_xyyyzz, g_y_0_zzzzz_xyyzzz, g_y_0_zzzzz_xyzzzz, g_y_0_zzzzz_xzzzzz, g_y_0_zzzzz_yyyyyy, g_y_0_zzzzz_yyyyyz, g_y_0_zzzzz_yyyyzz, g_y_0_zzzzz_yyyzzz, g_y_0_zzzzz_yyzzzz, g_y_0_zzzzz_yzzzzz, g_y_0_zzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzz_xxxxxx[k] = -g_y_0_zzzz_xxxxxx[k] * cd_z[k] + g_y_0_zzzz_xxxxxxz[k];

                g_y_0_zzzzz_xxxxxy[k] = -g_y_0_zzzz_xxxxxy[k] * cd_z[k] + g_y_0_zzzz_xxxxxyz[k];

                g_y_0_zzzzz_xxxxxz[k] = -g_y_0_zzzz_xxxxxz[k] * cd_z[k] + g_y_0_zzzz_xxxxxzz[k];

                g_y_0_zzzzz_xxxxyy[k] = -g_y_0_zzzz_xxxxyy[k] * cd_z[k] + g_y_0_zzzz_xxxxyyz[k];

                g_y_0_zzzzz_xxxxyz[k] = -g_y_0_zzzz_xxxxyz[k] * cd_z[k] + g_y_0_zzzz_xxxxyzz[k];

                g_y_0_zzzzz_xxxxzz[k] = -g_y_0_zzzz_xxxxzz[k] * cd_z[k] + g_y_0_zzzz_xxxxzzz[k];

                g_y_0_zzzzz_xxxyyy[k] = -g_y_0_zzzz_xxxyyy[k] * cd_z[k] + g_y_0_zzzz_xxxyyyz[k];

                g_y_0_zzzzz_xxxyyz[k] = -g_y_0_zzzz_xxxyyz[k] * cd_z[k] + g_y_0_zzzz_xxxyyzz[k];

                g_y_0_zzzzz_xxxyzz[k] = -g_y_0_zzzz_xxxyzz[k] * cd_z[k] + g_y_0_zzzz_xxxyzzz[k];

                g_y_0_zzzzz_xxxzzz[k] = -g_y_0_zzzz_xxxzzz[k] * cd_z[k] + g_y_0_zzzz_xxxzzzz[k];

                g_y_0_zzzzz_xxyyyy[k] = -g_y_0_zzzz_xxyyyy[k] * cd_z[k] + g_y_0_zzzz_xxyyyyz[k];

                g_y_0_zzzzz_xxyyyz[k] = -g_y_0_zzzz_xxyyyz[k] * cd_z[k] + g_y_0_zzzz_xxyyyzz[k];

                g_y_0_zzzzz_xxyyzz[k] = -g_y_0_zzzz_xxyyzz[k] * cd_z[k] + g_y_0_zzzz_xxyyzzz[k];

                g_y_0_zzzzz_xxyzzz[k] = -g_y_0_zzzz_xxyzzz[k] * cd_z[k] + g_y_0_zzzz_xxyzzzz[k];

                g_y_0_zzzzz_xxzzzz[k] = -g_y_0_zzzz_xxzzzz[k] * cd_z[k] + g_y_0_zzzz_xxzzzzz[k];

                g_y_0_zzzzz_xyyyyy[k] = -g_y_0_zzzz_xyyyyy[k] * cd_z[k] + g_y_0_zzzz_xyyyyyz[k];

                g_y_0_zzzzz_xyyyyz[k] = -g_y_0_zzzz_xyyyyz[k] * cd_z[k] + g_y_0_zzzz_xyyyyzz[k];

                g_y_0_zzzzz_xyyyzz[k] = -g_y_0_zzzz_xyyyzz[k] * cd_z[k] + g_y_0_zzzz_xyyyzzz[k];

                g_y_0_zzzzz_xyyzzz[k] = -g_y_0_zzzz_xyyzzz[k] * cd_z[k] + g_y_0_zzzz_xyyzzzz[k];

                g_y_0_zzzzz_xyzzzz[k] = -g_y_0_zzzz_xyzzzz[k] * cd_z[k] + g_y_0_zzzz_xyzzzzz[k];

                g_y_0_zzzzz_xzzzzz[k] = -g_y_0_zzzz_xzzzzz[k] * cd_z[k] + g_y_0_zzzz_xzzzzzz[k];

                g_y_0_zzzzz_yyyyyy[k] = -g_y_0_zzzz_yyyyyy[k] * cd_z[k] + g_y_0_zzzz_yyyyyyz[k];

                g_y_0_zzzzz_yyyyyz[k] = -g_y_0_zzzz_yyyyyz[k] * cd_z[k] + g_y_0_zzzz_yyyyyzz[k];

                g_y_0_zzzzz_yyyyzz[k] = -g_y_0_zzzz_yyyyzz[k] * cd_z[k] + g_y_0_zzzz_yyyyzzz[k];

                g_y_0_zzzzz_yyyzzz[k] = -g_y_0_zzzz_yyyzzz[k] * cd_z[k] + g_y_0_zzzz_yyyzzzz[k];

                g_y_0_zzzzz_yyzzzz[k] = -g_y_0_zzzz_yyzzzz[k] * cd_z[k] + g_y_0_zzzz_yyzzzzz[k];

                g_y_0_zzzzz_yzzzzz[k] = -g_y_0_zzzz_yzzzzz[k] * cd_z[k] + g_y_0_zzzz_yzzzzzz[k];

                g_y_0_zzzzz_zzzzzz[k] = -g_y_0_zzzz_zzzzzz[k] * cd_z[k] + g_y_0_zzzz_zzzzzzz[k];
            }
            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxx_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 0);

            auto g_z_0_xxxxx_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 1);

            auto g_z_0_xxxxx_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 2);

            auto g_z_0_xxxxx_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 3);

            auto g_z_0_xxxxx_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 4);

            auto g_z_0_xxxxx_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 5);

            auto g_z_0_xxxxx_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 6);

            auto g_z_0_xxxxx_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 7);

            auto g_z_0_xxxxx_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 8);

            auto g_z_0_xxxxx_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 9);

            auto g_z_0_xxxxx_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 10);

            auto g_z_0_xxxxx_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 11);

            auto g_z_0_xxxxx_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 12);

            auto g_z_0_xxxxx_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 13);

            auto g_z_0_xxxxx_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 14);

            auto g_z_0_xxxxx_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 15);

            auto g_z_0_xxxxx_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 16);

            auto g_z_0_xxxxx_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 17);

            auto g_z_0_xxxxx_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 18);

            auto g_z_0_xxxxx_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 19);

            auto g_z_0_xxxxx_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 20);

            auto g_z_0_xxxxx_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 21);

            auto g_z_0_xxxxx_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 22);

            auto g_z_0_xxxxx_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 23);

            auto g_z_0_xxxxx_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 24);

            auto g_z_0_xxxxx_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 25);

            auto g_z_0_xxxxx_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 26);

            auto g_z_0_xxxxx_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 27);

            #pragma omp simd aligned(cd_x, g_z_0_xxxx_xxxxxx, g_z_0_xxxx_xxxxxxx, g_z_0_xxxx_xxxxxxy, g_z_0_xxxx_xxxxxxz, g_z_0_xxxx_xxxxxy, g_z_0_xxxx_xxxxxyy, g_z_0_xxxx_xxxxxyz, g_z_0_xxxx_xxxxxz, g_z_0_xxxx_xxxxxzz, g_z_0_xxxx_xxxxyy, g_z_0_xxxx_xxxxyyy, g_z_0_xxxx_xxxxyyz, g_z_0_xxxx_xxxxyz, g_z_0_xxxx_xxxxyzz, g_z_0_xxxx_xxxxzz, g_z_0_xxxx_xxxxzzz, g_z_0_xxxx_xxxyyy, g_z_0_xxxx_xxxyyyy, g_z_0_xxxx_xxxyyyz, g_z_0_xxxx_xxxyyz, g_z_0_xxxx_xxxyyzz, g_z_0_xxxx_xxxyzz, g_z_0_xxxx_xxxyzzz, g_z_0_xxxx_xxxzzz, g_z_0_xxxx_xxxzzzz, g_z_0_xxxx_xxyyyy, g_z_0_xxxx_xxyyyyy, g_z_0_xxxx_xxyyyyz, g_z_0_xxxx_xxyyyz, g_z_0_xxxx_xxyyyzz, g_z_0_xxxx_xxyyzz, g_z_0_xxxx_xxyyzzz, g_z_0_xxxx_xxyzzz, g_z_0_xxxx_xxyzzzz, g_z_0_xxxx_xxzzzz, g_z_0_xxxx_xxzzzzz, g_z_0_xxxx_xyyyyy, g_z_0_xxxx_xyyyyyy, g_z_0_xxxx_xyyyyyz, g_z_0_xxxx_xyyyyz, g_z_0_xxxx_xyyyyzz, g_z_0_xxxx_xyyyzz, g_z_0_xxxx_xyyyzzz, g_z_0_xxxx_xyyzzz, g_z_0_xxxx_xyyzzzz, g_z_0_xxxx_xyzzzz, g_z_0_xxxx_xyzzzzz, g_z_0_xxxx_xzzzzz, g_z_0_xxxx_xzzzzzz, g_z_0_xxxx_yyyyyy, g_z_0_xxxx_yyyyyz, g_z_0_xxxx_yyyyzz, g_z_0_xxxx_yyyzzz, g_z_0_xxxx_yyzzzz, g_z_0_xxxx_yzzzzz, g_z_0_xxxx_zzzzzz, g_z_0_xxxxx_xxxxxx, g_z_0_xxxxx_xxxxxy, g_z_0_xxxxx_xxxxxz, g_z_0_xxxxx_xxxxyy, g_z_0_xxxxx_xxxxyz, g_z_0_xxxxx_xxxxzz, g_z_0_xxxxx_xxxyyy, g_z_0_xxxxx_xxxyyz, g_z_0_xxxxx_xxxyzz, g_z_0_xxxxx_xxxzzz, g_z_0_xxxxx_xxyyyy, g_z_0_xxxxx_xxyyyz, g_z_0_xxxxx_xxyyzz, g_z_0_xxxxx_xxyzzz, g_z_0_xxxxx_xxzzzz, g_z_0_xxxxx_xyyyyy, g_z_0_xxxxx_xyyyyz, g_z_0_xxxxx_xyyyzz, g_z_0_xxxxx_xyyzzz, g_z_0_xxxxx_xyzzzz, g_z_0_xxxxx_xzzzzz, g_z_0_xxxxx_yyyyyy, g_z_0_xxxxx_yyyyyz, g_z_0_xxxxx_yyyyzz, g_z_0_xxxxx_yyyzzz, g_z_0_xxxxx_yyzzzz, g_z_0_xxxxx_yzzzzz, g_z_0_xxxxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxx_xxxxxx[k] = -g_z_0_xxxx_xxxxxx[k] * cd_x[k] + g_z_0_xxxx_xxxxxxx[k];

                g_z_0_xxxxx_xxxxxy[k] = -g_z_0_xxxx_xxxxxy[k] * cd_x[k] + g_z_0_xxxx_xxxxxxy[k];

                g_z_0_xxxxx_xxxxxz[k] = -g_z_0_xxxx_xxxxxz[k] * cd_x[k] + g_z_0_xxxx_xxxxxxz[k];

                g_z_0_xxxxx_xxxxyy[k] = -g_z_0_xxxx_xxxxyy[k] * cd_x[k] + g_z_0_xxxx_xxxxxyy[k];

                g_z_0_xxxxx_xxxxyz[k] = -g_z_0_xxxx_xxxxyz[k] * cd_x[k] + g_z_0_xxxx_xxxxxyz[k];

                g_z_0_xxxxx_xxxxzz[k] = -g_z_0_xxxx_xxxxzz[k] * cd_x[k] + g_z_0_xxxx_xxxxxzz[k];

                g_z_0_xxxxx_xxxyyy[k] = -g_z_0_xxxx_xxxyyy[k] * cd_x[k] + g_z_0_xxxx_xxxxyyy[k];

                g_z_0_xxxxx_xxxyyz[k] = -g_z_0_xxxx_xxxyyz[k] * cd_x[k] + g_z_0_xxxx_xxxxyyz[k];

                g_z_0_xxxxx_xxxyzz[k] = -g_z_0_xxxx_xxxyzz[k] * cd_x[k] + g_z_0_xxxx_xxxxyzz[k];

                g_z_0_xxxxx_xxxzzz[k] = -g_z_0_xxxx_xxxzzz[k] * cd_x[k] + g_z_0_xxxx_xxxxzzz[k];

                g_z_0_xxxxx_xxyyyy[k] = -g_z_0_xxxx_xxyyyy[k] * cd_x[k] + g_z_0_xxxx_xxxyyyy[k];

                g_z_0_xxxxx_xxyyyz[k] = -g_z_0_xxxx_xxyyyz[k] * cd_x[k] + g_z_0_xxxx_xxxyyyz[k];

                g_z_0_xxxxx_xxyyzz[k] = -g_z_0_xxxx_xxyyzz[k] * cd_x[k] + g_z_0_xxxx_xxxyyzz[k];

                g_z_0_xxxxx_xxyzzz[k] = -g_z_0_xxxx_xxyzzz[k] * cd_x[k] + g_z_0_xxxx_xxxyzzz[k];

                g_z_0_xxxxx_xxzzzz[k] = -g_z_0_xxxx_xxzzzz[k] * cd_x[k] + g_z_0_xxxx_xxxzzzz[k];

                g_z_0_xxxxx_xyyyyy[k] = -g_z_0_xxxx_xyyyyy[k] * cd_x[k] + g_z_0_xxxx_xxyyyyy[k];

                g_z_0_xxxxx_xyyyyz[k] = -g_z_0_xxxx_xyyyyz[k] * cd_x[k] + g_z_0_xxxx_xxyyyyz[k];

                g_z_0_xxxxx_xyyyzz[k] = -g_z_0_xxxx_xyyyzz[k] * cd_x[k] + g_z_0_xxxx_xxyyyzz[k];

                g_z_0_xxxxx_xyyzzz[k] = -g_z_0_xxxx_xyyzzz[k] * cd_x[k] + g_z_0_xxxx_xxyyzzz[k];

                g_z_0_xxxxx_xyzzzz[k] = -g_z_0_xxxx_xyzzzz[k] * cd_x[k] + g_z_0_xxxx_xxyzzzz[k];

                g_z_0_xxxxx_xzzzzz[k] = -g_z_0_xxxx_xzzzzz[k] * cd_x[k] + g_z_0_xxxx_xxzzzzz[k];

                g_z_0_xxxxx_yyyyyy[k] = -g_z_0_xxxx_yyyyyy[k] * cd_x[k] + g_z_0_xxxx_xyyyyyy[k];

                g_z_0_xxxxx_yyyyyz[k] = -g_z_0_xxxx_yyyyyz[k] * cd_x[k] + g_z_0_xxxx_xyyyyyz[k];

                g_z_0_xxxxx_yyyyzz[k] = -g_z_0_xxxx_yyyyzz[k] * cd_x[k] + g_z_0_xxxx_xyyyyzz[k];

                g_z_0_xxxxx_yyyzzz[k] = -g_z_0_xxxx_yyyzzz[k] * cd_x[k] + g_z_0_xxxx_xyyyzzz[k];

                g_z_0_xxxxx_yyzzzz[k] = -g_z_0_xxxx_yyzzzz[k] * cd_x[k] + g_z_0_xxxx_xyyzzzz[k];

                g_z_0_xxxxx_yzzzzz[k] = -g_z_0_xxxx_yzzzzz[k] * cd_x[k] + g_z_0_xxxx_xyzzzzz[k];

                g_z_0_xxxxx_zzzzzz[k] = -g_z_0_xxxx_zzzzzz[k] * cd_x[k] + g_z_0_xxxx_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxy_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 28);

            auto g_z_0_xxxxy_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 29);

            auto g_z_0_xxxxy_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 30);

            auto g_z_0_xxxxy_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 31);

            auto g_z_0_xxxxy_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 32);

            auto g_z_0_xxxxy_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 33);

            auto g_z_0_xxxxy_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 34);

            auto g_z_0_xxxxy_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 35);

            auto g_z_0_xxxxy_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 36);

            auto g_z_0_xxxxy_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 37);

            auto g_z_0_xxxxy_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 38);

            auto g_z_0_xxxxy_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 39);

            auto g_z_0_xxxxy_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 40);

            auto g_z_0_xxxxy_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 41);

            auto g_z_0_xxxxy_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 42);

            auto g_z_0_xxxxy_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 43);

            auto g_z_0_xxxxy_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 44);

            auto g_z_0_xxxxy_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 45);

            auto g_z_0_xxxxy_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 46);

            auto g_z_0_xxxxy_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 47);

            auto g_z_0_xxxxy_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 48);

            auto g_z_0_xxxxy_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 49);

            auto g_z_0_xxxxy_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 50);

            auto g_z_0_xxxxy_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 51);

            auto g_z_0_xxxxy_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 52);

            auto g_z_0_xxxxy_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 53);

            auto g_z_0_xxxxy_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 54);

            auto g_z_0_xxxxy_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 55);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxy_xxxxxx, g_z_0_xxxxy_xxxxxy, g_z_0_xxxxy_xxxxxz, g_z_0_xxxxy_xxxxyy, g_z_0_xxxxy_xxxxyz, g_z_0_xxxxy_xxxxzz, g_z_0_xxxxy_xxxyyy, g_z_0_xxxxy_xxxyyz, g_z_0_xxxxy_xxxyzz, g_z_0_xxxxy_xxxzzz, g_z_0_xxxxy_xxyyyy, g_z_0_xxxxy_xxyyyz, g_z_0_xxxxy_xxyyzz, g_z_0_xxxxy_xxyzzz, g_z_0_xxxxy_xxzzzz, g_z_0_xxxxy_xyyyyy, g_z_0_xxxxy_xyyyyz, g_z_0_xxxxy_xyyyzz, g_z_0_xxxxy_xyyzzz, g_z_0_xxxxy_xyzzzz, g_z_0_xxxxy_xzzzzz, g_z_0_xxxxy_yyyyyy, g_z_0_xxxxy_yyyyyz, g_z_0_xxxxy_yyyyzz, g_z_0_xxxxy_yyyzzz, g_z_0_xxxxy_yyzzzz, g_z_0_xxxxy_yzzzzz, g_z_0_xxxxy_zzzzzz, g_z_0_xxxy_xxxxxx, g_z_0_xxxy_xxxxxxx, g_z_0_xxxy_xxxxxxy, g_z_0_xxxy_xxxxxxz, g_z_0_xxxy_xxxxxy, g_z_0_xxxy_xxxxxyy, g_z_0_xxxy_xxxxxyz, g_z_0_xxxy_xxxxxz, g_z_0_xxxy_xxxxxzz, g_z_0_xxxy_xxxxyy, g_z_0_xxxy_xxxxyyy, g_z_0_xxxy_xxxxyyz, g_z_0_xxxy_xxxxyz, g_z_0_xxxy_xxxxyzz, g_z_0_xxxy_xxxxzz, g_z_0_xxxy_xxxxzzz, g_z_0_xxxy_xxxyyy, g_z_0_xxxy_xxxyyyy, g_z_0_xxxy_xxxyyyz, g_z_0_xxxy_xxxyyz, g_z_0_xxxy_xxxyyzz, g_z_0_xxxy_xxxyzz, g_z_0_xxxy_xxxyzzz, g_z_0_xxxy_xxxzzz, g_z_0_xxxy_xxxzzzz, g_z_0_xxxy_xxyyyy, g_z_0_xxxy_xxyyyyy, g_z_0_xxxy_xxyyyyz, g_z_0_xxxy_xxyyyz, g_z_0_xxxy_xxyyyzz, g_z_0_xxxy_xxyyzz, g_z_0_xxxy_xxyyzzz, g_z_0_xxxy_xxyzzz, g_z_0_xxxy_xxyzzzz, g_z_0_xxxy_xxzzzz, g_z_0_xxxy_xxzzzzz, g_z_0_xxxy_xyyyyy, g_z_0_xxxy_xyyyyyy, g_z_0_xxxy_xyyyyyz, g_z_0_xxxy_xyyyyz, g_z_0_xxxy_xyyyyzz, g_z_0_xxxy_xyyyzz, g_z_0_xxxy_xyyyzzz, g_z_0_xxxy_xyyzzz, g_z_0_xxxy_xyyzzzz, g_z_0_xxxy_xyzzzz, g_z_0_xxxy_xyzzzzz, g_z_0_xxxy_xzzzzz, g_z_0_xxxy_xzzzzzz, g_z_0_xxxy_yyyyyy, g_z_0_xxxy_yyyyyz, g_z_0_xxxy_yyyyzz, g_z_0_xxxy_yyyzzz, g_z_0_xxxy_yyzzzz, g_z_0_xxxy_yzzzzz, g_z_0_xxxy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxy_xxxxxx[k] = -g_z_0_xxxy_xxxxxx[k] * cd_x[k] + g_z_0_xxxy_xxxxxxx[k];

                g_z_0_xxxxy_xxxxxy[k] = -g_z_0_xxxy_xxxxxy[k] * cd_x[k] + g_z_0_xxxy_xxxxxxy[k];

                g_z_0_xxxxy_xxxxxz[k] = -g_z_0_xxxy_xxxxxz[k] * cd_x[k] + g_z_0_xxxy_xxxxxxz[k];

                g_z_0_xxxxy_xxxxyy[k] = -g_z_0_xxxy_xxxxyy[k] * cd_x[k] + g_z_0_xxxy_xxxxxyy[k];

                g_z_0_xxxxy_xxxxyz[k] = -g_z_0_xxxy_xxxxyz[k] * cd_x[k] + g_z_0_xxxy_xxxxxyz[k];

                g_z_0_xxxxy_xxxxzz[k] = -g_z_0_xxxy_xxxxzz[k] * cd_x[k] + g_z_0_xxxy_xxxxxzz[k];

                g_z_0_xxxxy_xxxyyy[k] = -g_z_0_xxxy_xxxyyy[k] * cd_x[k] + g_z_0_xxxy_xxxxyyy[k];

                g_z_0_xxxxy_xxxyyz[k] = -g_z_0_xxxy_xxxyyz[k] * cd_x[k] + g_z_0_xxxy_xxxxyyz[k];

                g_z_0_xxxxy_xxxyzz[k] = -g_z_0_xxxy_xxxyzz[k] * cd_x[k] + g_z_0_xxxy_xxxxyzz[k];

                g_z_0_xxxxy_xxxzzz[k] = -g_z_0_xxxy_xxxzzz[k] * cd_x[k] + g_z_0_xxxy_xxxxzzz[k];

                g_z_0_xxxxy_xxyyyy[k] = -g_z_0_xxxy_xxyyyy[k] * cd_x[k] + g_z_0_xxxy_xxxyyyy[k];

                g_z_0_xxxxy_xxyyyz[k] = -g_z_0_xxxy_xxyyyz[k] * cd_x[k] + g_z_0_xxxy_xxxyyyz[k];

                g_z_0_xxxxy_xxyyzz[k] = -g_z_0_xxxy_xxyyzz[k] * cd_x[k] + g_z_0_xxxy_xxxyyzz[k];

                g_z_0_xxxxy_xxyzzz[k] = -g_z_0_xxxy_xxyzzz[k] * cd_x[k] + g_z_0_xxxy_xxxyzzz[k];

                g_z_0_xxxxy_xxzzzz[k] = -g_z_0_xxxy_xxzzzz[k] * cd_x[k] + g_z_0_xxxy_xxxzzzz[k];

                g_z_0_xxxxy_xyyyyy[k] = -g_z_0_xxxy_xyyyyy[k] * cd_x[k] + g_z_0_xxxy_xxyyyyy[k];

                g_z_0_xxxxy_xyyyyz[k] = -g_z_0_xxxy_xyyyyz[k] * cd_x[k] + g_z_0_xxxy_xxyyyyz[k];

                g_z_0_xxxxy_xyyyzz[k] = -g_z_0_xxxy_xyyyzz[k] * cd_x[k] + g_z_0_xxxy_xxyyyzz[k];

                g_z_0_xxxxy_xyyzzz[k] = -g_z_0_xxxy_xyyzzz[k] * cd_x[k] + g_z_0_xxxy_xxyyzzz[k];

                g_z_0_xxxxy_xyzzzz[k] = -g_z_0_xxxy_xyzzzz[k] * cd_x[k] + g_z_0_xxxy_xxyzzzz[k];

                g_z_0_xxxxy_xzzzzz[k] = -g_z_0_xxxy_xzzzzz[k] * cd_x[k] + g_z_0_xxxy_xxzzzzz[k];

                g_z_0_xxxxy_yyyyyy[k] = -g_z_0_xxxy_yyyyyy[k] * cd_x[k] + g_z_0_xxxy_xyyyyyy[k];

                g_z_0_xxxxy_yyyyyz[k] = -g_z_0_xxxy_yyyyyz[k] * cd_x[k] + g_z_0_xxxy_xyyyyyz[k];

                g_z_0_xxxxy_yyyyzz[k] = -g_z_0_xxxy_yyyyzz[k] * cd_x[k] + g_z_0_xxxy_xyyyyzz[k];

                g_z_0_xxxxy_yyyzzz[k] = -g_z_0_xxxy_yyyzzz[k] * cd_x[k] + g_z_0_xxxy_xyyyzzz[k];

                g_z_0_xxxxy_yyzzzz[k] = -g_z_0_xxxy_yyzzzz[k] * cd_x[k] + g_z_0_xxxy_xyyzzzz[k];

                g_z_0_xxxxy_yzzzzz[k] = -g_z_0_xxxy_yzzzzz[k] * cd_x[k] + g_z_0_xxxy_xyzzzzz[k];

                g_z_0_xxxxy_zzzzzz[k] = -g_z_0_xxxy_zzzzzz[k] * cd_x[k] + g_z_0_xxxy_xzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 56);

            auto g_z_0_xxxxz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 57);

            auto g_z_0_xxxxz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 58);

            auto g_z_0_xxxxz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 59);

            auto g_z_0_xxxxz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 60);

            auto g_z_0_xxxxz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 61);

            auto g_z_0_xxxxz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 62);

            auto g_z_0_xxxxz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 63);

            auto g_z_0_xxxxz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 64);

            auto g_z_0_xxxxz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 65);

            auto g_z_0_xxxxz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 66);

            auto g_z_0_xxxxz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 67);

            auto g_z_0_xxxxz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 68);

            auto g_z_0_xxxxz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 69);

            auto g_z_0_xxxxz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 70);

            auto g_z_0_xxxxz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 71);

            auto g_z_0_xxxxz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 72);

            auto g_z_0_xxxxz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 73);

            auto g_z_0_xxxxz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 74);

            auto g_z_0_xxxxz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 75);

            auto g_z_0_xxxxz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 76);

            auto g_z_0_xxxxz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 77);

            auto g_z_0_xxxxz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 78);

            auto g_z_0_xxxxz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 79);

            auto g_z_0_xxxxz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 80);

            auto g_z_0_xxxxz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 81);

            auto g_z_0_xxxxz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 82);

            auto g_z_0_xxxxz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxz_xxxxxx, g_z_0_xxxxz_xxxxxy, g_z_0_xxxxz_xxxxxz, g_z_0_xxxxz_xxxxyy, g_z_0_xxxxz_xxxxyz, g_z_0_xxxxz_xxxxzz, g_z_0_xxxxz_xxxyyy, g_z_0_xxxxz_xxxyyz, g_z_0_xxxxz_xxxyzz, g_z_0_xxxxz_xxxzzz, g_z_0_xxxxz_xxyyyy, g_z_0_xxxxz_xxyyyz, g_z_0_xxxxz_xxyyzz, g_z_0_xxxxz_xxyzzz, g_z_0_xxxxz_xxzzzz, g_z_0_xxxxz_xyyyyy, g_z_0_xxxxz_xyyyyz, g_z_0_xxxxz_xyyyzz, g_z_0_xxxxz_xyyzzz, g_z_0_xxxxz_xyzzzz, g_z_0_xxxxz_xzzzzz, g_z_0_xxxxz_yyyyyy, g_z_0_xxxxz_yyyyyz, g_z_0_xxxxz_yyyyzz, g_z_0_xxxxz_yyyzzz, g_z_0_xxxxz_yyzzzz, g_z_0_xxxxz_yzzzzz, g_z_0_xxxxz_zzzzzz, g_z_0_xxxz_xxxxxx, g_z_0_xxxz_xxxxxxx, g_z_0_xxxz_xxxxxxy, g_z_0_xxxz_xxxxxxz, g_z_0_xxxz_xxxxxy, g_z_0_xxxz_xxxxxyy, g_z_0_xxxz_xxxxxyz, g_z_0_xxxz_xxxxxz, g_z_0_xxxz_xxxxxzz, g_z_0_xxxz_xxxxyy, g_z_0_xxxz_xxxxyyy, g_z_0_xxxz_xxxxyyz, g_z_0_xxxz_xxxxyz, g_z_0_xxxz_xxxxyzz, g_z_0_xxxz_xxxxzz, g_z_0_xxxz_xxxxzzz, g_z_0_xxxz_xxxyyy, g_z_0_xxxz_xxxyyyy, g_z_0_xxxz_xxxyyyz, g_z_0_xxxz_xxxyyz, g_z_0_xxxz_xxxyyzz, g_z_0_xxxz_xxxyzz, g_z_0_xxxz_xxxyzzz, g_z_0_xxxz_xxxzzz, g_z_0_xxxz_xxxzzzz, g_z_0_xxxz_xxyyyy, g_z_0_xxxz_xxyyyyy, g_z_0_xxxz_xxyyyyz, g_z_0_xxxz_xxyyyz, g_z_0_xxxz_xxyyyzz, g_z_0_xxxz_xxyyzz, g_z_0_xxxz_xxyyzzz, g_z_0_xxxz_xxyzzz, g_z_0_xxxz_xxyzzzz, g_z_0_xxxz_xxzzzz, g_z_0_xxxz_xxzzzzz, g_z_0_xxxz_xyyyyy, g_z_0_xxxz_xyyyyyy, g_z_0_xxxz_xyyyyyz, g_z_0_xxxz_xyyyyz, g_z_0_xxxz_xyyyyzz, g_z_0_xxxz_xyyyzz, g_z_0_xxxz_xyyyzzz, g_z_0_xxxz_xyyzzz, g_z_0_xxxz_xyyzzzz, g_z_0_xxxz_xyzzzz, g_z_0_xxxz_xyzzzzz, g_z_0_xxxz_xzzzzz, g_z_0_xxxz_xzzzzzz, g_z_0_xxxz_yyyyyy, g_z_0_xxxz_yyyyyz, g_z_0_xxxz_yyyyzz, g_z_0_xxxz_yyyzzz, g_z_0_xxxz_yyzzzz, g_z_0_xxxz_yzzzzz, g_z_0_xxxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxz_xxxxxx[k] = -g_z_0_xxxz_xxxxxx[k] * cd_x[k] + g_z_0_xxxz_xxxxxxx[k];

                g_z_0_xxxxz_xxxxxy[k] = -g_z_0_xxxz_xxxxxy[k] * cd_x[k] + g_z_0_xxxz_xxxxxxy[k];

                g_z_0_xxxxz_xxxxxz[k] = -g_z_0_xxxz_xxxxxz[k] * cd_x[k] + g_z_0_xxxz_xxxxxxz[k];

                g_z_0_xxxxz_xxxxyy[k] = -g_z_0_xxxz_xxxxyy[k] * cd_x[k] + g_z_0_xxxz_xxxxxyy[k];

                g_z_0_xxxxz_xxxxyz[k] = -g_z_0_xxxz_xxxxyz[k] * cd_x[k] + g_z_0_xxxz_xxxxxyz[k];

                g_z_0_xxxxz_xxxxzz[k] = -g_z_0_xxxz_xxxxzz[k] * cd_x[k] + g_z_0_xxxz_xxxxxzz[k];

                g_z_0_xxxxz_xxxyyy[k] = -g_z_0_xxxz_xxxyyy[k] * cd_x[k] + g_z_0_xxxz_xxxxyyy[k];

                g_z_0_xxxxz_xxxyyz[k] = -g_z_0_xxxz_xxxyyz[k] * cd_x[k] + g_z_0_xxxz_xxxxyyz[k];

                g_z_0_xxxxz_xxxyzz[k] = -g_z_0_xxxz_xxxyzz[k] * cd_x[k] + g_z_0_xxxz_xxxxyzz[k];

                g_z_0_xxxxz_xxxzzz[k] = -g_z_0_xxxz_xxxzzz[k] * cd_x[k] + g_z_0_xxxz_xxxxzzz[k];

                g_z_0_xxxxz_xxyyyy[k] = -g_z_0_xxxz_xxyyyy[k] * cd_x[k] + g_z_0_xxxz_xxxyyyy[k];

                g_z_0_xxxxz_xxyyyz[k] = -g_z_0_xxxz_xxyyyz[k] * cd_x[k] + g_z_0_xxxz_xxxyyyz[k];

                g_z_0_xxxxz_xxyyzz[k] = -g_z_0_xxxz_xxyyzz[k] * cd_x[k] + g_z_0_xxxz_xxxyyzz[k];

                g_z_0_xxxxz_xxyzzz[k] = -g_z_0_xxxz_xxyzzz[k] * cd_x[k] + g_z_0_xxxz_xxxyzzz[k];

                g_z_0_xxxxz_xxzzzz[k] = -g_z_0_xxxz_xxzzzz[k] * cd_x[k] + g_z_0_xxxz_xxxzzzz[k];

                g_z_0_xxxxz_xyyyyy[k] = -g_z_0_xxxz_xyyyyy[k] * cd_x[k] + g_z_0_xxxz_xxyyyyy[k];

                g_z_0_xxxxz_xyyyyz[k] = -g_z_0_xxxz_xyyyyz[k] * cd_x[k] + g_z_0_xxxz_xxyyyyz[k];

                g_z_0_xxxxz_xyyyzz[k] = -g_z_0_xxxz_xyyyzz[k] * cd_x[k] + g_z_0_xxxz_xxyyyzz[k];

                g_z_0_xxxxz_xyyzzz[k] = -g_z_0_xxxz_xyyzzz[k] * cd_x[k] + g_z_0_xxxz_xxyyzzz[k];

                g_z_0_xxxxz_xyzzzz[k] = -g_z_0_xxxz_xyzzzz[k] * cd_x[k] + g_z_0_xxxz_xxyzzzz[k];

                g_z_0_xxxxz_xzzzzz[k] = -g_z_0_xxxz_xzzzzz[k] * cd_x[k] + g_z_0_xxxz_xxzzzzz[k];

                g_z_0_xxxxz_yyyyyy[k] = -g_z_0_xxxz_yyyyyy[k] * cd_x[k] + g_z_0_xxxz_xyyyyyy[k];

                g_z_0_xxxxz_yyyyyz[k] = -g_z_0_xxxz_yyyyyz[k] * cd_x[k] + g_z_0_xxxz_xyyyyyz[k];

                g_z_0_xxxxz_yyyyzz[k] = -g_z_0_xxxz_yyyyzz[k] * cd_x[k] + g_z_0_xxxz_xyyyyzz[k];

                g_z_0_xxxxz_yyyzzz[k] = -g_z_0_xxxz_yyyzzz[k] * cd_x[k] + g_z_0_xxxz_xyyyzzz[k];

                g_z_0_xxxxz_yyzzzz[k] = -g_z_0_xxxz_yyzzzz[k] * cd_x[k] + g_z_0_xxxz_xyyzzzz[k];

                g_z_0_xxxxz_yzzzzz[k] = -g_z_0_xxxz_yzzzzz[k] * cd_x[k] + g_z_0_xxxz_xyzzzzz[k];

                g_z_0_xxxxz_zzzzzz[k] = -g_z_0_xxxz_zzzzzz[k] * cd_x[k] + g_z_0_xxxz_xzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 84);

            auto g_z_0_xxxyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 85);

            auto g_z_0_xxxyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 86);

            auto g_z_0_xxxyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 87);

            auto g_z_0_xxxyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 88);

            auto g_z_0_xxxyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 89);

            auto g_z_0_xxxyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 90);

            auto g_z_0_xxxyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 91);

            auto g_z_0_xxxyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 92);

            auto g_z_0_xxxyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 93);

            auto g_z_0_xxxyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 94);

            auto g_z_0_xxxyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 95);

            auto g_z_0_xxxyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 96);

            auto g_z_0_xxxyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 97);

            auto g_z_0_xxxyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 98);

            auto g_z_0_xxxyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 99);

            auto g_z_0_xxxyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 100);

            auto g_z_0_xxxyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 101);

            auto g_z_0_xxxyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 102);

            auto g_z_0_xxxyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 103);

            auto g_z_0_xxxyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 104);

            auto g_z_0_xxxyy_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 105);

            auto g_z_0_xxxyy_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 106);

            auto g_z_0_xxxyy_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 107);

            auto g_z_0_xxxyy_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 108);

            auto g_z_0_xxxyy_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 109);

            auto g_z_0_xxxyy_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 110);

            auto g_z_0_xxxyy_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 111);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyy_xxxxxx, g_z_0_xxxyy_xxxxxy, g_z_0_xxxyy_xxxxxz, g_z_0_xxxyy_xxxxyy, g_z_0_xxxyy_xxxxyz, g_z_0_xxxyy_xxxxzz, g_z_0_xxxyy_xxxyyy, g_z_0_xxxyy_xxxyyz, g_z_0_xxxyy_xxxyzz, g_z_0_xxxyy_xxxzzz, g_z_0_xxxyy_xxyyyy, g_z_0_xxxyy_xxyyyz, g_z_0_xxxyy_xxyyzz, g_z_0_xxxyy_xxyzzz, g_z_0_xxxyy_xxzzzz, g_z_0_xxxyy_xyyyyy, g_z_0_xxxyy_xyyyyz, g_z_0_xxxyy_xyyyzz, g_z_0_xxxyy_xyyzzz, g_z_0_xxxyy_xyzzzz, g_z_0_xxxyy_xzzzzz, g_z_0_xxxyy_yyyyyy, g_z_0_xxxyy_yyyyyz, g_z_0_xxxyy_yyyyzz, g_z_0_xxxyy_yyyzzz, g_z_0_xxxyy_yyzzzz, g_z_0_xxxyy_yzzzzz, g_z_0_xxxyy_zzzzzz, g_z_0_xxyy_xxxxxx, g_z_0_xxyy_xxxxxxx, g_z_0_xxyy_xxxxxxy, g_z_0_xxyy_xxxxxxz, g_z_0_xxyy_xxxxxy, g_z_0_xxyy_xxxxxyy, g_z_0_xxyy_xxxxxyz, g_z_0_xxyy_xxxxxz, g_z_0_xxyy_xxxxxzz, g_z_0_xxyy_xxxxyy, g_z_0_xxyy_xxxxyyy, g_z_0_xxyy_xxxxyyz, g_z_0_xxyy_xxxxyz, g_z_0_xxyy_xxxxyzz, g_z_0_xxyy_xxxxzz, g_z_0_xxyy_xxxxzzz, g_z_0_xxyy_xxxyyy, g_z_0_xxyy_xxxyyyy, g_z_0_xxyy_xxxyyyz, g_z_0_xxyy_xxxyyz, g_z_0_xxyy_xxxyyzz, g_z_0_xxyy_xxxyzz, g_z_0_xxyy_xxxyzzz, g_z_0_xxyy_xxxzzz, g_z_0_xxyy_xxxzzzz, g_z_0_xxyy_xxyyyy, g_z_0_xxyy_xxyyyyy, g_z_0_xxyy_xxyyyyz, g_z_0_xxyy_xxyyyz, g_z_0_xxyy_xxyyyzz, g_z_0_xxyy_xxyyzz, g_z_0_xxyy_xxyyzzz, g_z_0_xxyy_xxyzzz, g_z_0_xxyy_xxyzzzz, g_z_0_xxyy_xxzzzz, g_z_0_xxyy_xxzzzzz, g_z_0_xxyy_xyyyyy, g_z_0_xxyy_xyyyyyy, g_z_0_xxyy_xyyyyyz, g_z_0_xxyy_xyyyyz, g_z_0_xxyy_xyyyyzz, g_z_0_xxyy_xyyyzz, g_z_0_xxyy_xyyyzzz, g_z_0_xxyy_xyyzzz, g_z_0_xxyy_xyyzzzz, g_z_0_xxyy_xyzzzz, g_z_0_xxyy_xyzzzzz, g_z_0_xxyy_xzzzzz, g_z_0_xxyy_xzzzzzz, g_z_0_xxyy_yyyyyy, g_z_0_xxyy_yyyyyz, g_z_0_xxyy_yyyyzz, g_z_0_xxyy_yyyzzz, g_z_0_xxyy_yyzzzz, g_z_0_xxyy_yzzzzz, g_z_0_xxyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyy_xxxxxx[k] = -g_z_0_xxyy_xxxxxx[k] * cd_x[k] + g_z_0_xxyy_xxxxxxx[k];

                g_z_0_xxxyy_xxxxxy[k] = -g_z_0_xxyy_xxxxxy[k] * cd_x[k] + g_z_0_xxyy_xxxxxxy[k];

                g_z_0_xxxyy_xxxxxz[k] = -g_z_0_xxyy_xxxxxz[k] * cd_x[k] + g_z_0_xxyy_xxxxxxz[k];

                g_z_0_xxxyy_xxxxyy[k] = -g_z_0_xxyy_xxxxyy[k] * cd_x[k] + g_z_0_xxyy_xxxxxyy[k];

                g_z_0_xxxyy_xxxxyz[k] = -g_z_0_xxyy_xxxxyz[k] * cd_x[k] + g_z_0_xxyy_xxxxxyz[k];

                g_z_0_xxxyy_xxxxzz[k] = -g_z_0_xxyy_xxxxzz[k] * cd_x[k] + g_z_0_xxyy_xxxxxzz[k];

                g_z_0_xxxyy_xxxyyy[k] = -g_z_0_xxyy_xxxyyy[k] * cd_x[k] + g_z_0_xxyy_xxxxyyy[k];

                g_z_0_xxxyy_xxxyyz[k] = -g_z_0_xxyy_xxxyyz[k] * cd_x[k] + g_z_0_xxyy_xxxxyyz[k];

                g_z_0_xxxyy_xxxyzz[k] = -g_z_0_xxyy_xxxyzz[k] * cd_x[k] + g_z_0_xxyy_xxxxyzz[k];

                g_z_0_xxxyy_xxxzzz[k] = -g_z_0_xxyy_xxxzzz[k] * cd_x[k] + g_z_0_xxyy_xxxxzzz[k];

                g_z_0_xxxyy_xxyyyy[k] = -g_z_0_xxyy_xxyyyy[k] * cd_x[k] + g_z_0_xxyy_xxxyyyy[k];

                g_z_0_xxxyy_xxyyyz[k] = -g_z_0_xxyy_xxyyyz[k] * cd_x[k] + g_z_0_xxyy_xxxyyyz[k];

                g_z_0_xxxyy_xxyyzz[k] = -g_z_0_xxyy_xxyyzz[k] * cd_x[k] + g_z_0_xxyy_xxxyyzz[k];

                g_z_0_xxxyy_xxyzzz[k] = -g_z_0_xxyy_xxyzzz[k] * cd_x[k] + g_z_0_xxyy_xxxyzzz[k];

                g_z_0_xxxyy_xxzzzz[k] = -g_z_0_xxyy_xxzzzz[k] * cd_x[k] + g_z_0_xxyy_xxxzzzz[k];

                g_z_0_xxxyy_xyyyyy[k] = -g_z_0_xxyy_xyyyyy[k] * cd_x[k] + g_z_0_xxyy_xxyyyyy[k];

                g_z_0_xxxyy_xyyyyz[k] = -g_z_0_xxyy_xyyyyz[k] * cd_x[k] + g_z_0_xxyy_xxyyyyz[k];

                g_z_0_xxxyy_xyyyzz[k] = -g_z_0_xxyy_xyyyzz[k] * cd_x[k] + g_z_0_xxyy_xxyyyzz[k];

                g_z_0_xxxyy_xyyzzz[k] = -g_z_0_xxyy_xyyzzz[k] * cd_x[k] + g_z_0_xxyy_xxyyzzz[k];

                g_z_0_xxxyy_xyzzzz[k] = -g_z_0_xxyy_xyzzzz[k] * cd_x[k] + g_z_0_xxyy_xxyzzzz[k];

                g_z_0_xxxyy_xzzzzz[k] = -g_z_0_xxyy_xzzzzz[k] * cd_x[k] + g_z_0_xxyy_xxzzzzz[k];

                g_z_0_xxxyy_yyyyyy[k] = -g_z_0_xxyy_yyyyyy[k] * cd_x[k] + g_z_0_xxyy_xyyyyyy[k];

                g_z_0_xxxyy_yyyyyz[k] = -g_z_0_xxyy_yyyyyz[k] * cd_x[k] + g_z_0_xxyy_xyyyyyz[k];

                g_z_0_xxxyy_yyyyzz[k] = -g_z_0_xxyy_yyyyzz[k] * cd_x[k] + g_z_0_xxyy_xyyyyzz[k];

                g_z_0_xxxyy_yyyzzz[k] = -g_z_0_xxyy_yyyzzz[k] * cd_x[k] + g_z_0_xxyy_xyyyzzz[k];

                g_z_0_xxxyy_yyzzzz[k] = -g_z_0_xxyy_yyzzzz[k] * cd_x[k] + g_z_0_xxyy_xyyzzzz[k];

                g_z_0_xxxyy_yzzzzz[k] = -g_z_0_xxyy_yzzzzz[k] * cd_x[k] + g_z_0_xxyy_xyzzzzz[k];

                g_z_0_xxxyy_zzzzzz[k] = -g_z_0_xxyy_zzzzzz[k] * cd_x[k] + g_z_0_xxyy_xzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 112);

            auto g_z_0_xxxyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 113);

            auto g_z_0_xxxyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 114);

            auto g_z_0_xxxyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 115);

            auto g_z_0_xxxyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 116);

            auto g_z_0_xxxyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 117);

            auto g_z_0_xxxyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 118);

            auto g_z_0_xxxyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 119);

            auto g_z_0_xxxyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 120);

            auto g_z_0_xxxyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 121);

            auto g_z_0_xxxyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 122);

            auto g_z_0_xxxyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 123);

            auto g_z_0_xxxyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 124);

            auto g_z_0_xxxyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 125);

            auto g_z_0_xxxyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 126);

            auto g_z_0_xxxyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 127);

            auto g_z_0_xxxyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 128);

            auto g_z_0_xxxyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 129);

            auto g_z_0_xxxyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 130);

            auto g_z_0_xxxyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 131);

            auto g_z_0_xxxyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 132);

            auto g_z_0_xxxyz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 133);

            auto g_z_0_xxxyz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 134);

            auto g_z_0_xxxyz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 135);

            auto g_z_0_xxxyz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 136);

            auto g_z_0_xxxyz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 137);

            auto g_z_0_xxxyz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 138);

            auto g_z_0_xxxyz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 139);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyz_xxxxxx, g_z_0_xxxyz_xxxxxy, g_z_0_xxxyz_xxxxxz, g_z_0_xxxyz_xxxxyy, g_z_0_xxxyz_xxxxyz, g_z_0_xxxyz_xxxxzz, g_z_0_xxxyz_xxxyyy, g_z_0_xxxyz_xxxyyz, g_z_0_xxxyz_xxxyzz, g_z_0_xxxyz_xxxzzz, g_z_0_xxxyz_xxyyyy, g_z_0_xxxyz_xxyyyz, g_z_0_xxxyz_xxyyzz, g_z_0_xxxyz_xxyzzz, g_z_0_xxxyz_xxzzzz, g_z_0_xxxyz_xyyyyy, g_z_0_xxxyz_xyyyyz, g_z_0_xxxyz_xyyyzz, g_z_0_xxxyz_xyyzzz, g_z_0_xxxyz_xyzzzz, g_z_0_xxxyz_xzzzzz, g_z_0_xxxyz_yyyyyy, g_z_0_xxxyz_yyyyyz, g_z_0_xxxyz_yyyyzz, g_z_0_xxxyz_yyyzzz, g_z_0_xxxyz_yyzzzz, g_z_0_xxxyz_yzzzzz, g_z_0_xxxyz_zzzzzz, g_z_0_xxyz_xxxxxx, g_z_0_xxyz_xxxxxxx, g_z_0_xxyz_xxxxxxy, g_z_0_xxyz_xxxxxxz, g_z_0_xxyz_xxxxxy, g_z_0_xxyz_xxxxxyy, g_z_0_xxyz_xxxxxyz, g_z_0_xxyz_xxxxxz, g_z_0_xxyz_xxxxxzz, g_z_0_xxyz_xxxxyy, g_z_0_xxyz_xxxxyyy, g_z_0_xxyz_xxxxyyz, g_z_0_xxyz_xxxxyz, g_z_0_xxyz_xxxxyzz, g_z_0_xxyz_xxxxzz, g_z_0_xxyz_xxxxzzz, g_z_0_xxyz_xxxyyy, g_z_0_xxyz_xxxyyyy, g_z_0_xxyz_xxxyyyz, g_z_0_xxyz_xxxyyz, g_z_0_xxyz_xxxyyzz, g_z_0_xxyz_xxxyzz, g_z_0_xxyz_xxxyzzz, g_z_0_xxyz_xxxzzz, g_z_0_xxyz_xxxzzzz, g_z_0_xxyz_xxyyyy, g_z_0_xxyz_xxyyyyy, g_z_0_xxyz_xxyyyyz, g_z_0_xxyz_xxyyyz, g_z_0_xxyz_xxyyyzz, g_z_0_xxyz_xxyyzz, g_z_0_xxyz_xxyyzzz, g_z_0_xxyz_xxyzzz, g_z_0_xxyz_xxyzzzz, g_z_0_xxyz_xxzzzz, g_z_0_xxyz_xxzzzzz, g_z_0_xxyz_xyyyyy, g_z_0_xxyz_xyyyyyy, g_z_0_xxyz_xyyyyyz, g_z_0_xxyz_xyyyyz, g_z_0_xxyz_xyyyyzz, g_z_0_xxyz_xyyyzz, g_z_0_xxyz_xyyyzzz, g_z_0_xxyz_xyyzzz, g_z_0_xxyz_xyyzzzz, g_z_0_xxyz_xyzzzz, g_z_0_xxyz_xyzzzzz, g_z_0_xxyz_xzzzzz, g_z_0_xxyz_xzzzzzz, g_z_0_xxyz_yyyyyy, g_z_0_xxyz_yyyyyz, g_z_0_xxyz_yyyyzz, g_z_0_xxyz_yyyzzz, g_z_0_xxyz_yyzzzz, g_z_0_xxyz_yzzzzz, g_z_0_xxyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyz_xxxxxx[k] = -g_z_0_xxyz_xxxxxx[k] * cd_x[k] + g_z_0_xxyz_xxxxxxx[k];

                g_z_0_xxxyz_xxxxxy[k] = -g_z_0_xxyz_xxxxxy[k] * cd_x[k] + g_z_0_xxyz_xxxxxxy[k];

                g_z_0_xxxyz_xxxxxz[k] = -g_z_0_xxyz_xxxxxz[k] * cd_x[k] + g_z_0_xxyz_xxxxxxz[k];

                g_z_0_xxxyz_xxxxyy[k] = -g_z_0_xxyz_xxxxyy[k] * cd_x[k] + g_z_0_xxyz_xxxxxyy[k];

                g_z_0_xxxyz_xxxxyz[k] = -g_z_0_xxyz_xxxxyz[k] * cd_x[k] + g_z_0_xxyz_xxxxxyz[k];

                g_z_0_xxxyz_xxxxzz[k] = -g_z_0_xxyz_xxxxzz[k] * cd_x[k] + g_z_0_xxyz_xxxxxzz[k];

                g_z_0_xxxyz_xxxyyy[k] = -g_z_0_xxyz_xxxyyy[k] * cd_x[k] + g_z_0_xxyz_xxxxyyy[k];

                g_z_0_xxxyz_xxxyyz[k] = -g_z_0_xxyz_xxxyyz[k] * cd_x[k] + g_z_0_xxyz_xxxxyyz[k];

                g_z_0_xxxyz_xxxyzz[k] = -g_z_0_xxyz_xxxyzz[k] * cd_x[k] + g_z_0_xxyz_xxxxyzz[k];

                g_z_0_xxxyz_xxxzzz[k] = -g_z_0_xxyz_xxxzzz[k] * cd_x[k] + g_z_0_xxyz_xxxxzzz[k];

                g_z_0_xxxyz_xxyyyy[k] = -g_z_0_xxyz_xxyyyy[k] * cd_x[k] + g_z_0_xxyz_xxxyyyy[k];

                g_z_0_xxxyz_xxyyyz[k] = -g_z_0_xxyz_xxyyyz[k] * cd_x[k] + g_z_0_xxyz_xxxyyyz[k];

                g_z_0_xxxyz_xxyyzz[k] = -g_z_0_xxyz_xxyyzz[k] * cd_x[k] + g_z_0_xxyz_xxxyyzz[k];

                g_z_0_xxxyz_xxyzzz[k] = -g_z_0_xxyz_xxyzzz[k] * cd_x[k] + g_z_0_xxyz_xxxyzzz[k];

                g_z_0_xxxyz_xxzzzz[k] = -g_z_0_xxyz_xxzzzz[k] * cd_x[k] + g_z_0_xxyz_xxxzzzz[k];

                g_z_0_xxxyz_xyyyyy[k] = -g_z_0_xxyz_xyyyyy[k] * cd_x[k] + g_z_0_xxyz_xxyyyyy[k];

                g_z_0_xxxyz_xyyyyz[k] = -g_z_0_xxyz_xyyyyz[k] * cd_x[k] + g_z_0_xxyz_xxyyyyz[k];

                g_z_0_xxxyz_xyyyzz[k] = -g_z_0_xxyz_xyyyzz[k] * cd_x[k] + g_z_0_xxyz_xxyyyzz[k];

                g_z_0_xxxyz_xyyzzz[k] = -g_z_0_xxyz_xyyzzz[k] * cd_x[k] + g_z_0_xxyz_xxyyzzz[k];

                g_z_0_xxxyz_xyzzzz[k] = -g_z_0_xxyz_xyzzzz[k] * cd_x[k] + g_z_0_xxyz_xxyzzzz[k];

                g_z_0_xxxyz_xzzzzz[k] = -g_z_0_xxyz_xzzzzz[k] * cd_x[k] + g_z_0_xxyz_xxzzzzz[k];

                g_z_0_xxxyz_yyyyyy[k] = -g_z_0_xxyz_yyyyyy[k] * cd_x[k] + g_z_0_xxyz_xyyyyyy[k];

                g_z_0_xxxyz_yyyyyz[k] = -g_z_0_xxyz_yyyyyz[k] * cd_x[k] + g_z_0_xxyz_xyyyyyz[k];

                g_z_0_xxxyz_yyyyzz[k] = -g_z_0_xxyz_yyyyzz[k] * cd_x[k] + g_z_0_xxyz_xyyyyzz[k];

                g_z_0_xxxyz_yyyzzz[k] = -g_z_0_xxyz_yyyzzz[k] * cd_x[k] + g_z_0_xxyz_xyyyzzz[k];

                g_z_0_xxxyz_yyzzzz[k] = -g_z_0_xxyz_yyzzzz[k] * cd_x[k] + g_z_0_xxyz_xyyzzzz[k];

                g_z_0_xxxyz_yzzzzz[k] = -g_z_0_xxyz_yzzzzz[k] * cd_x[k] + g_z_0_xxyz_xyzzzzz[k];

                g_z_0_xxxyz_zzzzzz[k] = -g_z_0_xxyz_zzzzzz[k] * cd_x[k] + g_z_0_xxyz_xzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 140);

            auto g_z_0_xxxzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 141);

            auto g_z_0_xxxzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 142);

            auto g_z_0_xxxzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 143);

            auto g_z_0_xxxzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 144);

            auto g_z_0_xxxzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 145);

            auto g_z_0_xxxzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 146);

            auto g_z_0_xxxzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 147);

            auto g_z_0_xxxzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 148);

            auto g_z_0_xxxzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 149);

            auto g_z_0_xxxzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 150);

            auto g_z_0_xxxzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 151);

            auto g_z_0_xxxzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 152);

            auto g_z_0_xxxzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 153);

            auto g_z_0_xxxzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 154);

            auto g_z_0_xxxzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 155);

            auto g_z_0_xxxzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 156);

            auto g_z_0_xxxzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 157);

            auto g_z_0_xxxzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 158);

            auto g_z_0_xxxzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 159);

            auto g_z_0_xxxzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 160);

            auto g_z_0_xxxzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 161);

            auto g_z_0_xxxzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 162);

            auto g_z_0_xxxzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 163);

            auto g_z_0_xxxzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 164);

            auto g_z_0_xxxzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 165);

            auto g_z_0_xxxzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 166);

            auto g_z_0_xxxzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_x, g_z_0_xxxzz_xxxxxx, g_z_0_xxxzz_xxxxxy, g_z_0_xxxzz_xxxxxz, g_z_0_xxxzz_xxxxyy, g_z_0_xxxzz_xxxxyz, g_z_0_xxxzz_xxxxzz, g_z_0_xxxzz_xxxyyy, g_z_0_xxxzz_xxxyyz, g_z_0_xxxzz_xxxyzz, g_z_0_xxxzz_xxxzzz, g_z_0_xxxzz_xxyyyy, g_z_0_xxxzz_xxyyyz, g_z_0_xxxzz_xxyyzz, g_z_0_xxxzz_xxyzzz, g_z_0_xxxzz_xxzzzz, g_z_0_xxxzz_xyyyyy, g_z_0_xxxzz_xyyyyz, g_z_0_xxxzz_xyyyzz, g_z_0_xxxzz_xyyzzz, g_z_0_xxxzz_xyzzzz, g_z_0_xxxzz_xzzzzz, g_z_0_xxxzz_yyyyyy, g_z_0_xxxzz_yyyyyz, g_z_0_xxxzz_yyyyzz, g_z_0_xxxzz_yyyzzz, g_z_0_xxxzz_yyzzzz, g_z_0_xxxzz_yzzzzz, g_z_0_xxxzz_zzzzzz, g_z_0_xxzz_xxxxxx, g_z_0_xxzz_xxxxxxx, g_z_0_xxzz_xxxxxxy, g_z_0_xxzz_xxxxxxz, g_z_0_xxzz_xxxxxy, g_z_0_xxzz_xxxxxyy, g_z_0_xxzz_xxxxxyz, g_z_0_xxzz_xxxxxz, g_z_0_xxzz_xxxxxzz, g_z_0_xxzz_xxxxyy, g_z_0_xxzz_xxxxyyy, g_z_0_xxzz_xxxxyyz, g_z_0_xxzz_xxxxyz, g_z_0_xxzz_xxxxyzz, g_z_0_xxzz_xxxxzz, g_z_0_xxzz_xxxxzzz, g_z_0_xxzz_xxxyyy, g_z_0_xxzz_xxxyyyy, g_z_0_xxzz_xxxyyyz, g_z_0_xxzz_xxxyyz, g_z_0_xxzz_xxxyyzz, g_z_0_xxzz_xxxyzz, g_z_0_xxzz_xxxyzzz, g_z_0_xxzz_xxxzzz, g_z_0_xxzz_xxxzzzz, g_z_0_xxzz_xxyyyy, g_z_0_xxzz_xxyyyyy, g_z_0_xxzz_xxyyyyz, g_z_0_xxzz_xxyyyz, g_z_0_xxzz_xxyyyzz, g_z_0_xxzz_xxyyzz, g_z_0_xxzz_xxyyzzz, g_z_0_xxzz_xxyzzz, g_z_0_xxzz_xxyzzzz, g_z_0_xxzz_xxzzzz, g_z_0_xxzz_xxzzzzz, g_z_0_xxzz_xyyyyy, g_z_0_xxzz_xyyyyyy, g_z_0_xxzz_xyyyyyz, g_z_0_xxzz_xyyyyz, g_z_0_xxzz_xyyyyzz, g_z_0_xxzz_xyyyzz, g_z_0_xxzz_xyyyzzz, g_z_0_xxzz_xyyzzz, g_z_0_xxzz_xyyzzzz, g_z_0_xxzz_xyzzzz, g_z_0_xxzz_xyzzzzz, g_z_0_xxzz_xzzzzz, g_z_0_xxzz_xzzzzzz, g_z_0_xxzz_yyyyyy, g_z_0_xxzz_yyyyyz, g_z_0_xxzz_yyyyzz, g_z_0_xxzz_yyyzzz, g_z_0_xxzz_yyzzzz, g_z_0_xxzz_yzzzzz, g_z_0_xxzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzz_xxxxxx[k] = -g_z_0_xxzz_xxxxxx[k] * cd_x[k] + g_z_0_xxzz_xxxxxxx[k];

                g_z_0_xxxzz_xxxxxy[k] = -g_z_0_xxzz_xxxxxy[k] * cd_x[k] + g_z_0_xxzz_xxxxxxy[k];

                g_z_0_xxxzz_xxxxxz[k] = -g_z_0_xxzz_xxxxxz[k] * cd_x[k] + g_z_0_xxzz_xxxxxxz[k];

                g_z_0_xxxzz_xxxxyy[k] = -g_z_0_xxzz_xxxxyy[k] * cd_x[k] + g_z_0_xxzz_xxxxxyy[k];

                g_z_0_xxxzz_xxxxyz[k] = -g_z_0_xxzz_xxxxyz[k] * cd_x[k] + g_z_0_xxzz_xxxxxyz[k];

                g_z_0_xxxzz_xxxxzz[k] = -g_z_0_xxzz_xxxxzz[k] * cd_x[k] + g_z_0_xxzz_xxxxxzz[k];

                g_z_0_xxxzz_xxxyyy[k] = -g_z_0_xxzz_xxxyyy[k] * cd_x[k] + g_z_0_xxzz_xxxxyyy[k];

                g_z_0_xxxzz_xxxyyz[k] = -g_z_0_xxzz_xxxyyz[k] * cd_x[k] + g_z_0_xxzz_xxxxyyz[k];

                g_z_0_xxxzz_xxxyzz[k] = -g_z_0_xxzz_xxxyzz[k] * cd_x[k] + g_z_0_xxzz_xxxxyzz[k];

                g_z_0_xxxzz_xxxzzz[k] = -g_z_0_xxzz_xxxzzz[k] * cd_x[k] + g_z_0_xxzz_xxxxzzz[k];

                g_z_0_xxxzz_xxyyyy[k] = -g_z_0_xxzz_xxyyyy[k] * cd_x[k] + g_z_0_xxzz_xxxyyyy[k];

                g_z_0_xxxzz_xxyyyz[k] = -g_z_0_xxzz_xxyyyz[k] * cd_x[k] + g_z_0_xxzz_xxxyyyz[k];

                g_z_0_xxxzz_xxyyzz[k] = -g_z_0_xxzz_xxyyzz[k] * cd_x[k] + g_z_0_xxzz_xxxyyzz[k];

                g_z_0_xxxzz_xxyzzz[k] = -g_z_0_xxzz_xxyzzz[k] * cd_x[k] + g_z_0_xxzz_xxxyzzz[k];

                g_z_0_xxxzz_xxzzzz[k] = -g_z_0_xxzz_xxzzzz[k] * cd_x[k] + g_z_0_xxzz_xxxzzzz[k];

                g_z_0_xxxzz_xyyyyy[k] = -g_z_0_xxzz_xyyyyy[k] * cd_x[k] + g_z_0_xxzz_xxyyyyy[k];

                g_z_0_xxxzz_xyyyyz[k] = -g_z_0_xxzz_xyyyyz[k] * cd_x[k] + g_z_0_xxzz_xxyyyyz[k];

                g_z_0_xxxzz_xyyyzz[k] = -g_z_0_xxzz_xyyyzz[k] * cd_x[k] + g_z_0_xxzz_xxyyyzz[k];

                g_z_0_xxxzz_xyyzzz[k] = -g_z_0_xxzz_xyyzzz[k] * cd_x[k] + g_z_0_xxzz_xxyyzzz[k];

                g_z_0_xxxzz_xyzzzz[k] = -g_z_0_xxzz_xyzzzz[k] * cd_x[k] + g_z_0_xxzz_xxyzzzz[k];

                g_z_0_xxxzz_xzzzzz[k] = -g_z_0_xxzz_xzzzzz[k] * cd_x[k] + g_z_0_xxzz_xxzzzzz[k];

                g_z_0_xxxzz_yyyyyy[k] = -g_z_0_xxzz_yyyyyy[k] * cd_x[k] + g_z_0_xxzz_xyyyyyy[k];

                g_z_0_xxxzz_yyyyyz[k] = -g_z_0_xxzz_yyyyyz[k] * cd_x[k] + g_z_0_xxzz_xyyyyyz[k];

                g_z_0_xxxzz_yyyyzz[k] = -g_z_0_xxzz_yyyyzz[k] * cd_x[k] + g_z_0_xxzz_xyyyyzz[k];

                g_z_0_xxxzz_yyyzzz[k] = -g_z_0_xxzz_yyyzzz[k] * cd_x[k] + g_z_0_xxzz_xyyyzzz[k];

                g_z_0_xxxzz_yyzzzz[k] = -g_z_0_xxzz_yyzzzz[k] * cd_x[k] + g_z_0_xxzz_xyyzzzz[k];

                g_z_0_xxxzz_yzzzzz[k] = -g_z_0_xxzz_yzzzzz[k] * cd_x[k] + g_z_0_xxzz_xyzzzzz[k];

                g_z_0_xxxzz_zzzzzz[k] = -g_z_0_xxzz_zzzzzz[k] * cd_x[k] + g_z_0_xxzz_xzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 168);

            auto g_z_0_xxyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 169);

            auto g_z_0_xxyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 170);

            auto g_z_0_xxyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 171);

            auto g_z_0_xxyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 172);

            auto g_z_0_xxyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 173);

            auto g_z_0_xxyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 174);

            auto g_z_0_xxyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 175);

            auto g_z_0_xxyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 176);

            auto g_z_0_xxyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 177);

            auto g_z_0_xxyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 178);

            auto g_z_0_xxyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 179);

            auto g_z_0_xxyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 180);

            auto g_z_0_xxyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 181);

            auto g_z_0_xxyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 182);

            auto g_z_0_xxyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 183);

            auto g_z_0_xxyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 184);

            auto g_z_0_xxyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 185);

            auto g_z_0_xxyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 186);

            auto g_z_0_xxyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 187);

            auto g_z_0_xxyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 188);

            auto g_z_0_xxyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 189);

            auto g_z_0_xxyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 190);

            auto g_z_0_xxyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 191);

            auto g_z_0_xxyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 192);

            auto g_z_0_xxyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 193);

            auto g_z_0_xxyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 194);

            auto g_z_0_xxyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 195);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyy_xxxxxx, g_z_0_xxyyy_xxxxxy, g_z_0_xxyyy_xxxxxz, g_z_0_xxyyy_xxxxyy, g_z_0_xxyyy_xxxxyz, g_z_0_xxyyy_xxxxzz, g_z_0_xxyyy_xxxyyy, g_z_0_xxyyy_xxxyyz, g_z_0_xxyyy_xxxyzz, g_z_0_xxyyy_xxxzzz, g_z_0_xxyyy_xxyyyy, g_z_0_xxyyy_xxyyyz, g_z_0_xxyyy_xxyyzz, g_z_0_xxyyy_xxyzzz, g_z_0_xxyyy_xxzzzz, g_z_0_xxyyy_xyyyyy, g_z_0_xxyyy_xyyyyz, g_z_0_xxyyy_xyyyzz, g_z_0_xxyyy_xyyzzz, g_z_0_xxyyy_xyzzzz, g_z_0_xxyyy_xzzzzz, g_z_0_xxyyy_yyyyyy, g_z_0_xxyyy_yyyyyz, g_z_0_xxyyy_yyyyzz, g_z_0_xxyyy_yyyzzz, g_z_0_xxyyy_yyzzzz, g_z_0_xxyyy_yzzzzz, g_z_0_xxyyy_zzzzzz, g_z_0_xyyy_xxxxxx, g_z_0_xyyy_xxxxxxx, g_z_0_xyyy_xxxxxxy, g_z_0_xyyy_xxxxxxz, g_z_0_xyyy_xxxxxy, g_z_0_xyyy_xxxxxyy, g_z_0_xyyy_xxxxxyz, g_z_0_xyyy_xxxxxz, g_z_0_xyyy_xxxxxzz, g_z_0_xyyy_xxxxyy, g_z_0_xyyy_xxxxyyy, g_z_0_xyyy_xxxxyyz, g_z_0_xyyy_xxxxyz, g_z_0_xyyy_xxxxyzz, g_z_0_xyyy_xxxxzz, g_z_0_xyyy_xxxxzzz, g_z_0_xyyy_xxxyyy, g_z_0_xyyy_xxxyyyy, g_z_0_xyyy_xxxyyyz, g_z_0_xyyy_xxxyyz, g_z_0_xyyy_xxxyyzz, g_z_0_xyyy_xxxyzz, g_z_0_xyyy_xxxyzzz, g_z_0_xyyy_xxxzzz, g_z_0_xyyy_xxxzzzz, g_z_0_xyyy_xxyyyy, g_z_0_xyyy_xxyyyyy, g_z_0_xyyy_xxyyyyz, g_z_0_xyyy_xxyyyz, g_z_0_xyyy_xxyyyzz, g_z_0_xyyy_xxyyzz, g_z_0_xyyy_xxyyzzz, g_z_0_xyyy_xxyzzz, g_z_0_xyyy_xxyzzzz, g_z_0_xyyy_xxzzzz, g_z_0_xyyy_xxzzzzz, g_z_0_xyyy_xyyyyy, g_z_0_xyyy_xyyyyyy, g_z_0_xyyy_xyyyyyz, g_z_0_xyyy_xyyyyz, g_z_0_xyyy_xyyyyzz, g_z_0_xyyy_xyyyzz, g_z_0_xyyy_xyyyzzz, g_z_0_xyyy_xyyzzz, g_z_0_xyyy_xyyzzzz, g_z_0_xyyy_xyzzzz, g_z_0_xyyy_xyzzzzz, g_z_0_xyyy_xzzzzz, g_z_0_xyyy_xzzzzzz, g_z_0_xyyy_yyyyyy, g_z_0_xyyy_yyyyyz, g_z_0_xyyy_yyyyzz, g_z_0_xyyy_yyyzzz, g_z_0_xyyy_yyzzzz, g_z_0_xyyy_yzzzzz, g_z_0_xyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyy_xxxxxx[k] = -g_z_0_xyyy_xxxxxx[k] * cd_x[k] + g_z_0_xyyy_xxxxxxx[k];

                g_z_0_xxyyy_xxxxxy[k] = -g_z_0_xyyy_xxxxxy[k] * cd_x[k] + g_z_0_xyyy_xxxxxxy[k];

                g_z_0_xxyyy_xxxxxz[k] = -g_z_0_xyyy_xxxxxz[k] * cd_x[k] + g_z_0_xyyy_xxxxxxz[k];

                g_z_0_xxyyy_xxxxyy[k] = -g_z_0_xyyy_xxxxyy[k] * cd_x[k] + g_z_0_xyyy_xxxxxyy[k];

                g_z_0_xxyyy_xxxxyz[k] = -g_z_0_xyyy_xxxxyz[k] * cd_x[k] + g_z_0_xyyy_xxxxxyz[k];

                g_z_0_xxyyy_xxxxzz[k] = -g_z_0_xyyy_xxxxzz[k] * cd_x[k] + g_z_0_xyyy_xxxxxzz[k];

                g_z_0_xxyyy_xxxyyy[k] = -g_z_0_xyyy_xxxyyy[k] * cd_x[k] + g_z_0_xyyy_xxxxyyy[k];

                g_z_0_xxyyy_xxxyyz[k] = -g_z_0_xyyy_xxxyyz[k] * cd_x[k] + g_z_0_xyyy_xxxxyyz[k];

                g_z_0_xxyyy_xxxyzz[k] = -g_z_0_xyyy_xxxyzz[k] * cd_x[k] + g_z_0_xyyy_xxxxyzz[k];

                g_z_0_xxyyy_xxxzzz[k] = -g_z_0_xyyy_xxxzzz[k] * cd_x[k] + g_z_0_xyyy_xxxxzzz[k];

                g_z_0_xxyyy_xxyyyy[k] = -g_z_0_xyyy_xxyyyy[k] * cd_x[k] + g_z_0_xyyy_xxxyyyy[k];

                g_z_0_xxyyy_xxyyyz[k] = -g_z_0_xyyy_xxyyyz[k] * cd_x[k] + g_z_0_xyyy_xxxyyyz[k];

                g_z_0_xxyyy_xxyyzz[k] = -g_z_0_xyyy_xxyyzz[k] * cd_x[k] + g_z_0_xyyy_xxxyyzz[k];

                g_z_0_xxyyy_xxyzzz[k] = -g_z_0_xyyy_xxyzzz[k] * cd_x[k] + g_z_0_xyyy_xxxyzzz[k];

                g_z_0_xxyyy_xxzzzz[k] = -g_z_0_xyyy_xxzzzz[k] * cd_x[k] + g_z_0_xyyy_xxxzzzz[k];

                g_z_0_xxyyy_xyyyyy[k] = -g_z_0_xyyy_xyyyyy[k] * cd_x[k] + g_z_0_xyyy_xxyyyyy[k];

                g_z_0_xxyyy_xyyyyz[k] = -g_z_0_xyyy_xyyyyz[k] * cd_x[k] + g_z_0_xyyy_xxyyyyz[k];

                g_z_0_xxyyy_xyyyzz[k] = -g_z_0_xyyy_xyyyzz[k] * cd_x[k] + g_z_0_xyyy_xxyyyzz[k];

                g_z_0_xxyyy_xyyzzz[k] = -g_z_0_xyyy_xyyzzz[k] * cd_x[k] + g_z_0_xyyy_xxyyzzz[k];

                g_z_0_xxyyy_xyzzzz[k] = -g_z_0_xyyy_xyzzzz[k] * cd_x[k] + g_z_0_xyyy_xxyzzzz[k];

                g_z_0_xxyyy_xzzzzz[k] = -g_z_0_xyyy_xzzzzz[k] * cd_x[k] + g_z_0_xyyy_xxzzzzz[k];

                g_z_0_xxyyy_yyyyyy[k] = -g_z_0_xyyy_yyyyyy[k] * cd_x[k] + g_z_0_xyyy_xyyyyyy[k];

                g_z_0_xxyyy_yyyyyz[k] = -g_z_0_xyyy_yyyyyz[k] * cd_x[k] + g_z_0_xyyy_xyyyyyz[k];

                g_z_0_xxyyy_yyyyzz[k] = -g_z_0_xyyy_yyyyzz[k] * cd_x[k] + g_z_0_xyyy_xyyyyzz[k];

                g_z_0_xxyyy_yyyzzz[k] = -g_z_0_xyyy_yyyzzz[k] * cd_x[k] + g_z_0_xyyy_xyyyzzz[k];

                g_z_0_xxyyy_yyzzzz[k] = -g_z_0_xyyy_yyzzzz[k] * cd_x[k] + g_z_0_xyyy_xyyzzzz[k];

                g_z_0_xxyyy_yzzzzz[k] = -g_z_0_xyyy_yzzzzz[k] * cd_x[k] + g_z_0_xyyy_xyzzzzz[k];

                g_z_0_xxyyy_zzzzzz[k] = -g_z_0_xyyy_zzzzzz[k] * cd_x[k] + g_z_0_xyyy_xzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 196);

            auto g_z_0_xxyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 197);

            auto g_z_0_xxyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 198);

            auto g_z_0_xxyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 199);

            auto g_z_0_xxyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 200);

            auto g_z_0_xxyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 201);

            auto g_z_0_xxyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 202);

            auto g_z_0_xxyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 203);

            auto g_z_0_xxyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 204);

            auto g_z_0_xxyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 205);

            auto g_z_0_xxyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 206);

            auto g_z_0_xxyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 207);

            auto g_z_0_xxyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 208);

            auto g_z_0_xxyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 209);

            auto g_z_0_xxyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 210);

            auto g_z_0_xxyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 211);

            auto g_z_0_xxyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 212);

            auto g_z_0_xxyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 213);

            auto g_z_0_xxyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 214);

            auto g_z_0_xxyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 215);

            auto g_z_0_xxyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 216);

            auto g_z_0_xxyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 217);

            auto g_z_0_xxyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 218);

            auto g_z_0_xxyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 219);

            auto g_z_0_xxyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 220);

            auto g_z_0_xxyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 221);

            auto g_z_0_xxyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 222);

            auto g_z_0_xxyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 223);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyz_xxxxxx, g_z_0_xxyyz_xxxxxy, g_z_0_xxyyz_xxxxxz, g_z_0_xxyyz_xxxxyy, g_z_0_xxyyz_xxxxyz, g_z_0_xxyyz_xxxxzz, g_z_0_xxyyz_xxxyyy, g_z_0_xxyyz_xxxyyz, g_z_0_xxyyz_xxxyzz, g_z_0_xxyyz_xxxzzz, g_z_0_xxyyz_xxyyyy, g_z_0_xxyyz_xxyyyz, g_z_0_xxyyz_xxyyzz, g_z_0_xxyyz_xxyzzz, g_z_0_xxyyz_xxzzzz, g_z_0_xxyyz_xyyyyy, g_z_0_xxyyz_xyyyyz, g_z_0_xxyyz_xyyyzz, g_z_0_xxyyz_xyyzzz, g_z_0_xxyyz_xyzzzz, g_z_0_xxyyz_xzzzzz, g_z_0_xxyyz_yyyyyy, g_z_0_xxyyz_yyyyyz, g_z_0_xxyyz_yyyyzz, g_z_0_xxyyz_yyyzzz, g_z_0_xxyyz_yyzzzz, g_z_0_xxyyz_yzzzzz, g_z_0_xxyyz_zzzzzz, g_z_0_xyyz_xxxxxx, g_z_0_xyyz_xxxxxxx, g_z_0_xyyz_xxxxxxy, g_z_0_xyyz_xxxxxxz, g_z_0_xyyz_xxxxxy, g_z_0_xyyz_xxxxxyy, g_z_0_xyyz_xxxxxyz, g_z_0_xyyz_xxxxxz, g_z_0_xyyz_xxxxxzz, g_z_0_xyyz_xxxxyy, g_z_0_xyyz_xxxxyyy, g_z_0_xyyz_xxxxyyz, g_z_0_xyyz_xxxxyz, g_z_0_xyyz_xxxxyzz, g_z_0_xyyz_xxxxzz, g_z_0_xyyz_xxxxzzz, g_z_0_xyyz_xxxyyy, g_z_0_xyyz_xxxyyyy, g_z_0_xyyz_xxxyyyz, g_z_0_xyyz_xxxyyz, g_z_0_xyyz_xxxyyzz, g_z_0_xyyz_xxxyzz, g_z_0_xyyz_xxxyzzz, g_z_0_xyyz_xxxzzz, g_z_0_xyyz_xxxzzzz, g_z_0_xyyz_xxyyyy, g_z_0_xyyz_xxyyyyy, g_z_0_xyyz_xxyyyyz, g_z_0_xyyz_xxyyyz, g_z_0_xyyz_xxyyyzz, g_z_0_xyyz_xxyyzz, g_z_0_xyyz_xxyyzzz, g_z_0_xyyz_xxyzzz, g_z_0_xyyz_xxyzzzz, g_z_0_xyyz_xxzzzz, g_z_0_xyyz_xxzzzzz, g_z_0_xyyz_xyyyyy, g_z_0_xyyz_xyyyyyy, g_z_0_xyyz_xyyyyyz, g_z_0_xyyz_xyyyyz, g_z_0_xyyz_xyyyyzz, g_z_0_xyyz_xyyyzz, g_z_0_xyyz_xyyyzzz, g_z_0_xyyz_xyyzzz, g_z_0_xyyz_xyyzzzz, g_z_0_xyyz_xyzzzz, g_z_0_xyyz_xyzzzzz, g_z_0_xyyz_xzzzzz, g_z_0_xyyz_xzzzzzz, g_z_0_xyyz_yyyyyy, g_z_0_xyyz_yyyyyz, g_z_0_xyyz_yyyyzz, g_z_0_xyyz_yyyzzz, g_z_0_xyyz_yyzzzz, g_z_0_xyyz_yzzzzz, g_z_0_xyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyz_xxxxxx[k] = -g_z_0_xyyz_xxxxxx[k] * cd_x[k] + g_z_0_xyyz_xxxxxxx[k];

                g_z_0_xxyyz_xxxxxy[k] = -g_z_0_xyyz_xxxxxy[k] * cd_x[k] + g_z_0_xyyz_xxxxxxy[k];

                g_z_0_xxyyz_xxxxxz[k] = -g_z_0_xyyz_xxxxxz[k] * cd_x[k] + g_z_0_xyyz_xxxxxxz[k];

                g_z_0_xxyyz_xxxxyy[k] = -g_z_0_xyyz_xxxxyy[k] * cd_x[k] + g_z_0_xyyz_xxxxxyy[k];

                g_z_0_xxyyz_xxxxyz[k] = -g_z_0_xyyz_xxxxyz[k] * cd_x[k] + g_z_0_xyyz_xxxxxyz[k];

                g_z_0_xxyyz_xxxxzz[k] = -g_z_0_xyyz_xxxxzz[k] * cd_x[k] + g_z_0_xyyz_xxxxxzz[k];

                g_z_0_xxyyz_xxxyyy[k] = -g_z_0_xyyz_xxxyyy[k] * cd_x[k] + g_z_0_xyyz_xxxxyyy[k];

                g_z_0_xxyyz_xxxyyz[k] = -g_z_0_xyyz_xxxyyz[k] * cd_x[k] + g_z_0_xyyz_xxxxyyz[k];

                g_z_0_xxyyz_xxxyzz[k] = -g_z_0_xyyz_xxxyzz[k] * cd_x[k] + g_z_0_xyyz_xxxxyzz[k];

                g_z_0_xxyyz_xxxzzz[k] = -g_z_0_xyyz_xxxzzz[k] * cd_x[k] + g_z_0_xyyz_xxxxzzz[k];

                g_z_0_xxyyz_xxyyyy[k] = -g_z_0_xyyz_xxyyyy[k] * cd_x[k] + g_z_0_xyyz_xxxyyyy[k];

                g_z_0_xxyyz_xxyyyz[k] = -g_z_0_xyyz_xxyyyz[k] * cd_x[k] + g_z_0_xyyz_xxxyyyz[k];

                g_z_0_xxyyz_xxyyzz[k] = -g_z_0_xyyz_xxyyzz[k] * cd_x[k] + g_z_0_xyyz_xxxyyzz[k];

                g_z_0_xxyyz_xxyzzz[k] = -g_z_0_xyyz_xxyzzz[k] * cd_x[k] + g_z_0_xyyz_xxxyzzz[k];

                g_z_0_xxyyz_xxzzzz[k] = -g_z_0_xyyz_xxzzzz[k] * cd_x[k] + g_z_0_xyyz_xxxzzzz[k];

                g_z_0_xxyyz_xyyyyy[k] = -g_z_0_xyyz_xyyyyy[k] * cd_x[k] + g_z_0_xyyz_xxyyyyy[k];

                g_z_0_xxyyz_xyyyyz[k] = -g_z_0_xyyz_xyyyyz[k] * cd_x[k] + g_z_0_xyyz_xxyyyyz[k];

                g_z_0_xxyyz_xyyyzz[k] = -g_z_0_xyyz_xyyyzz[k] * cd_x[k] + g_z_0_xyyz_xxyyyzz[k];

                g_z_0_xxyyz_xyyzzz[k] = -g_z_0_xyyz_xyyzzz[k] * cd_x[k] + g_z_0_xyyz_xxyyzzz[k];

                g_z_0_xxyyz_xyzzzz[k] = -g_z_0_xyyz_xyzzzz[k] * cd_x[k] + g_z_0_xyyz_xxyzzzz[k];

                g_z_0_xxyyz_xzzzzz[k] = -g_z_0_xyyz_xzzzzz[k] * cd_x[k] + g_z_0_xyyz_xxzzzzz[k];

                g_z_0_xxyyz_yyyyyy[k] = -g_z_0_xyyz_yyyyyy[k] * cd_x[k] + g_z_0_xyyz_xyyyyyy[k];

                g_z_0_xxyyz_yyyyyz[k] = -g_z_0_xyyz_yyyyyz[k] * cd_x[k] + g_z_0_xyyz_xyyyyyz[k];

                g_z_0_xxyyz_yyyyzz[k] = -g_z_0_xyyz_yyyyzz[k] * cd_x[k] + g_z_0_xyyz_xyyyyzz[k];

                g_z_0_xxyyz_yyyzzz[k] = -g_z_0_xyyz_yyyzzz[k] * cd_x[k] + g_z_0_xyyz_xyyyzzz[k];

                g_z_0_xxyyz_yyzzzz[k] = -g_z_0_xyyz_yyzzzz[k] * cd_x[k] + g_z_0_xyyz_xyyzzzz[k];

                g_z_0_xxyyz_yzzzzz[k] = -g_z_0_xyyz_yzzzzz[k] * cd_x[k] + g_z_0_xyyz_xyzzzzz[k];

                g_z_0_xxyyz_zzzzzz[k] = -g_z_0_xyyz_zzzzzz[k] * cd_x[k] + g_z_0_xyyz_xzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 224);

            auto g_z_0_xxyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 225);

            auto g_z_0_xxyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 226);

            auto g_z_0_xxyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 227);

            auto g_z_0_xxyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 228);

            auto g_z_0_xxyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 229);

            auto g_z_0_xxyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 230);

            auto g_z_0_xxyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 231);

            auto g_z_0_xxyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 232);

            auto g_z_0_xxyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 233);

            auto g_z_0_xxyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 234);

            auto g_z_0_xxyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 235);

            auto g_z_0_xxyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 236);

            auto g_z_0_xxyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 237);

            auto g_z_0_xxyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 238);

            auto g_z_0_xxyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 239);

            auto g_z_0_xxyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 240);

            auto g_z_0_xxyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 241);

            auto g_z_0_xxyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 242);

            auto g_z_0_xxyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 243);

            auto g_z_0_xxyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 244);

            auto g_z_0_xxyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 245);

            auto g_z_0_xxyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 246);

            auto g_z_0_xxyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 247);

            auto g_z_0_xxyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 248);

            auto g_z_0_xxyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 249);

            auto g_z_0_xxyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 250);

            auto g_z_0_xxyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 251);

            #pragma omp simd aligned(cd_x, g_z_0_xxyzz_xxxxxx, g_z_0_xxyzz_xxxxxy, g_z_0_xxyzz_xxxxxz, g_z_0_xxyzz_xxxxyy, g_z_0_xxyzz_xxxxyz, g_z_0_xxyzz_xxxxzz, g_z_0_xxyzz_xxxyyy, g_z_0_xxyzz_xxxyyz, g_z_0_xxyzz_xxxyzz, g_z_0_xxyzz_xxxzzz, g_z_0_xxyzz_xxyyyy, g_z_0_xxyzz_xxyyyz, g_z_0_xxyzz_xxyyzz, g_z_0_xxyzz_xxyzzz, g_z_0_xxyzz_xxzzzz, g_z_0_xxyzz_xyyyyy, g_z_0_xxyzz_xyyyyz, g_z_0_xxyzz_xyyyzz, g_z_0_xxyzz_xyyzzz, g_z_0_xxyzz_xyzzzz, g_z_0_xxyzz_xzzzzz, g_z_0_xxyzz_yyyyyy, g_z_0_xxyzz_yyyyyz, g_z_0_xxyzz_yyyyzz, g_z_0_xxyzz_yyyzzz, g_z_0_xxyzz_yyzzzz, g_z_0_xxyzz_yzzzzz, g_z_0_xxyzz_zzzzzz, g_z_0_xyzz_xxxxxx, g_z_0_xyzz_xxxxxxx, g_z_0_xyzz_xxxxxxy, g_z_0_xyzz_xxxxxxz, g_z_0_xyzz_xxxxxy, g_z_0_xyzz_xxxxxyy, g_z_0_xyzz_xxxxxyz, g_z_0_xyzz_xxxxxz, g_z_0_xyzz_xxxxxzz, g_z_0_xyzz_xxxxyy, g_z_0_xyzz_xxxxyyy, g_z_0_xyzz_xxxxyyz, g_z_0_xyzz_xxxxyz, g_z_0_xyzz_xxxxyzz, g_z_0_xyzz_xxxxzz, g_z_0_xyzz_xxxxzzz, g_z_0_xyzz_xxxyyy, g_z_0_xyzz_xxxyyyy, g_z_0_xyzz_xxxyyyz, g_z_0_xyzz_xxxyyz, g_z_0_xyzz_xxxyyzz, g_z_0_xyzz_xxxyzz, g_z_0_xyzz_xxxyzzz, g_z_0_xyzz_xxxzzz, g_z_0_xyzz_xxxzzzz, g_z_0_xyzz_xxyyyy, g_z_0_xyzz_xxyyyyy, g_z_0_xyzz_xxyyyyz, g_z_0_xyzz_xxyyyz, g_z_0_xyzz_xxyyyzz, g_z_0_xyzz_xxyyzz, g_z_0_xyzz_xxyyzzz, g_z_0_xyzz_xxyzzz, g_z_0_xyzz_xxyzzzz, g_z_0_xyzz_xxzzzz, g_z_0_xyzz_xxzzzzz, g_z_0_xyzz_xyyyyy, g_z_0_xyzz_xyyyyyy, g_z_0_xyzz_xyyyyyz, g_z_0_xyzz_xyyyyz, g_z_0_xyzz_xyyyyzz, g_z_0_xyzz_xyyyzz, g_z_0_xyzz_xyyyzzz, g_z_0_xyzz_xyyzzz, g_z_0_xyzz_xyyzzzz, g_z_0_xyzz_xyzzzz, g_z_0_xyzz_xyzzzzz, g_z_0_xyzz_xzzzzz, g_z_0_xyzz_xzzzzzz, g_z_0_xyzz_yyyyyy, g_z_0_xyzz_yyyyyz, g_z_0_xyzz_yyyyzz, g_z_0_xyzz_yyyzzz, g_z_0_xyzz_yyzzzz, g_z_0_xyzz_yzzzzz, g_z_0_xyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzz_xxxxxx[k] = -g_z_0_xyzz_xxxxxx[k] * cd_x[k] + g_z_0_xyzz_xxxxxxx[k];

                g_z_0_xxyzz_xxxxxy[k] = -g_z_0_xyzz_xxxxxy[k] * cd_x[k] + g_z_0_xyzz_xxxxxxy[k];

                g_z_0_xxyzz_xxxxxz[k] = -g_z_0_xyzz_xxxxxz[k] * cd_x[k] + g_z_0_xyzz_xxxxxxz[k];

                g_z_0_xxyzz_xxxxyy[k] = -g_z_0_xyzz_xxxxyy[k] * cd_x[k] + g_z_0_xyzz_xxxxxyy[k];

                g_z_0_xxyzz_xxxxyz[k] = -g_z_0_xyzz_xxxxyz[k] * cd_x[k] + g_z_0_xyzz_xxxxxyz[k];

                g_z_0_xxyzz_xxxxzz[k] = -g_z_0_xyzz_xxxxzz[k] * cd_x[k] + g_z_0_xyzz_xxxxxzz[k];

                g_z_0_xxyzz_xxxyyy[k] = -g_z_0_xyzz_xxxyyy[k] * cd_x[k] + g_z_0_xyzz_xxxxyyy[k];

                g_z_0_xxyzz_xxxyyz[k] = -g_z_0_xyzz_xxxyyz[k] * cd_x[k] + g_z_0_xyzz_xxxxyyz[k];

                g_z_0_xxyzz_xxxyzz[k] = -g_z_0_xyzz_xxxyzz[k] * cd_x[k] + g_z_0_xyzz_xxxxyzz[k];

                g_z_0_xxyzz_xxxzzz[k] = -g_z_0_xyzz_xxxzzz[k] * cd_x[k] + g_z_0_xyzz_xxxxzzz[k];

                g_z_0_xxyzz_xxyyyy[k] = -g_z_0_xyzz_xxyyyy[k] * cd_x[k] + g_z_0_xyzz_xxxyyyy[k];

                g_z_0_xxyzz_xxyyyz[k] = -g_z_0_xyzz_xxyyyz[k] * cd_x[k] + g_z_0_xyzz_xxxyyyz[k];

                g_z_0_xxyzz_xxyyzz[k] = -g_z_0_xyzz_xxyyzz[k] * cd_x[k] + g_z_0_xyzz_xxxyyzz[k];

                g_z_0_xxyzz_xxyzzz[k] = -g_z_0_xyzz_xxyzzz[k] * cd_x[k] + g_z_0_xyzz_xxxyzzz[k];

                g_z_0_xxyzz_xxzzzz[k] = -g_z_0_xyzz_xxzzzz[k] * cd_x[k] + g_z_0_xyzz_xxxzzzz[k];

                g_z_0_xxyzz_xyyyyy[k] = -g_z_0_xyzz_xyyyyy[k] * cd_x[k] + g_z_0_xyzz_xxyyyyy[k];

                g_z_0_xxyzz_xyyyyz[k] = -g_z_0_xyzz_xyyyyz[k] * cd_x[k] + g_z_0_xyzz_xxyyyyz[k];

                g_z_0_xxyzz_xyyyzz[k] = -g_z_0_xyzz_xyyyzz[k] * cd_x[k] + g_z_0_xyzz_xxyyyzz[k];

                g_z_0_xxyzz_xyyzzz[k] = -g_z_0_xyzz_xyyzzz[k] * cd_x[k] + g_z_0_xyzz_xxyyzzz[k];

                g_z_0_xxyzz_xyzzzz[k] = -g_z_0_xyzz_xyzzzz[k] * cd_x[k] + g_z_0_xyzz_xxyzzzz[k];

                g_z_0_xxyzz_xzzzzz[k] = -g_z_0_xyzz_xzzzzz[k] * cd_x[k] + g_z_0_xyzz_xxzzzzz[k];

                g_z_0_xxyzz_yyyyyy[k] = -g_z_0_xyzz_yyyyyy[k] * cd_x[k] + g_z_0_xyzz_xyyyyyy[k];

                g_z_0_xxyzz_yyyyyz[k] = -g_z_0_xyzz_yyyyyz[k] * cd_x[k] + g_z_0_xyzz_xyyyyyz[k];

                g_z_0_xxyzz_yyyyzz[k] = -g_z_0_xyzz_yyyyzz[k] * cd_x[k] + g_z_0_xyzz_xyyyyzz[k];

                g_z_0_xxyzz_yyyzzz[k] = -g_z_0_xyzz_yyyzzz[k] * cd_x[k] + g_z_0_xyzz_xyyyzzz[k];

                g_z_0_xxyzz_yyzzzz[k] = -g_z_0_xyzz_yyzzzz[k] * cd_x[k] + g_z_0_xyzz_xyyzzzz[k];

                g_z_0_xxyzz_yzzzzz[k] = -g_z_0_xyzz_yzzzzz[k] * cd_x[k] + g_z_0_xyzz_xyzzzzz[k];

                g_z_0_xxyzz_zzzzzz[k] = -g_z_0_xyzz_zzzzzz[k] * cd_x[k] + g_z_0_xyzz_xzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 252);

            auto g_z_0_xxzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 253);

            auto g_z_0_xxzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 254);

            auto g_z_0_xxzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 255);

            auto g_z_0_xxzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 256);

            auto g_z_0_xxzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 257);

            auto g_z_0_xxzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 258);

            auto g_z_0_xxzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 259);

            auto g_z_0_xxzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 260);

            auto g_z_0_xxzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 261);

            auto g_z_0_xxzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 262);

            auto g_z_0_xxzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 263);

            auto g_z_0_xxzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 264);

            auto g_z_0_xxzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 265);

            auto g_z_0_xxzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 266);

            auto g_z_0_xxzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 267);

            auto g_z_0_xxzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 268);

            auto g_z_0_xxzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 269);

            auto g_z_0_xxzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 270);

            auto g_z_0_xxzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 271);

            auto g_z_0_xxzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 272);

            auto g_z_0_xxzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 273);

            auto g_z_0_xxzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 274);

            auto g_z_0_xxzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 275);

            auto g_z_0_xxzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 276);

            auto g_z_0_xxzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 277);

            auto g_z_0_xxzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 278);

            auto g_z_0_xxzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 279);

            #pragma omp simd aligned(cd_x, g_z_0_xxzzz_xxxxxx, g_z_0_xxzzz_xxxxxy, g_z_0_xxzzz_xxxxxz, g_z_0_xxzzz_xxxxyy, g_z_0_xxzzz_xxxxyz, g_z_0_xxzzz_xxxxzz, g_z_0_xxzzz_xxxyyy, g_z_0_xxzzz_xxxyyz, g_z_0_xxzzz_xxxyzz, g_z_0_xxzzz_xxxzzz, g_z_0_xxzzz_xxyyyy, g_z_0_xxzzz_xxyyyz, g_z_0_xxzzz_xxyyzz, g_z_0_xxzzz_xxyzzz, g_z_0_xxzzz_xxzzzz, g_z_0_xxzzz_xyyyyy, g_z_0_xxzzz_xyyyyz, g_z_0_xxzzz_xyyyzz, g_z_0_xxzzz_xyyzzz, g_z_0_xxzzz_xyzzzz, g_z_0_xxzzz_xzzzzz, g_z_0_xxzzz_yyyyyy, g_z_0_xxzzz_yyyyyz, g_z_0_xxzzz_yyyyzz, g_z_0_xxzzz_yyyzzz, g_z_0_xxzzz_yyzzzz, g_z_0_xxzzz_yzzzzz, g_z_0_xxzzz_zzzzzz, g_z_0_xzzz_xxxxxx, g_z_0_xzzz_xxxxxxx, g_z_0_xzzz_xxxxxxy, g_z_0_xzzz_xxxxxxz, g_z_0_xzzz_xxxxxy, g_z_0_xzzz_xxxxxyy, g_z_0_xzzz_xxxxxyz, g_z_0_xzzz_xxxxxz, g_z_0_xzzz_xxxxxzz, g_z_0_xzzz_xxxxyy, g_z_0_xzzz_xxxxyyy, g_z_0_xzzz_xxxxyyz, g_z_0_xzzz_xxxxyz, g_z_0_xzzz_xxxxyzz, g_z_0_xzzz_xxxxzz, g_z_0_xzzz_xxxxzzz, g_z_0_xzzz_xxxyyy, g_z_0_xzzz_xxxyyyy, g_z_0_xzzz_xxxyyyz, g_z_0_xzzz_xxxyyz, g_z_0_xzzz_xxxyyzz, g_z_0_xzzz_xxxyzz, g_z_0_xzzz_xxxyzzz, g_z_0_xzzz_xxxzzz, g_z_0_xzzz_xxxzzzz, g_z_0_xzzz_xxyyyy, g_z_0_xzzz_xxyyyyy, g_z_0_xzzz_xxyyyyz, g_z_0_xzzz_xxyyyz, g_z_0_xzzz_xxyyyzz, g_z_0_xzzz_xxyyzz, g_z_0_xzzz_xxyyzzz, g_z_0_xzzz_xxyzzz, g_z_0_xzzz_xxyzzzz, g_z_0_xzzz_xxzzzz, g_z_0_xzzz_xxzzzzz, g_z_0_xzzz_xyyyyy, g_z_0_xzzz_xyyyyyy, g_z_0_xzzz_xyyyyyz, g_z_0_xzzz_xyyyyz, g_z_0_xzzz_xyyyyzz, g_z_0_xzzz_xyyyzz, g_z_0_xzzz_xyyyzzz, g_z_0_xzzz_xyyzzz, g_z_0_xzzz_xyyzzzz, g_z_0_xzzz_xyzzzz, g_z_0_xzzz_xyzzzzz, g_z_0_xzzz_xzzzzz, g_z_0_xzzz_xzzzzzz, g_z_0_xzzz_yyyyyy, g_z_0_xzzz_yyyyyz, g_z_0_xzzz_yyyyzz, g_z_0_xzzz_yyyzzz, g_z_0_xzzz_yyzzzz, g_z_0_xzzz_yzzzzz, g_z_0_xzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzz_xxxxxx[k] = -g_z_0_xzzz_xxxxxx[k] * cd_x[k] + g_z_0_xzzz_xxxxxxx[k];

                g_z_0_xxzzz_xxxxxy[k] = -g_z_0_xzzz_xxxxxy[k] * cd_x[k] + g_z_0_xzzz_xxxxxxy[k];

                g_z_0_xxzzz_xxxxxz[k] = -g_z_0_xzzz_xxxxxz[k] * cd_x[k] + g_z_0_xzzz_xxxxxxz[k];

                g_z_0_xxzzz_xxxxyy[k] = -g_z_0_xzzz_xxxxyy[k] * cd_x[k] + g_z_0_xzzz_xxxxxyy[k];

                g_z_0_xxzzz_xxxxyz[k] = -g_z_0_xzzz_xxxxyz[k] * cd_x[k] + g_z_0_xzzz_xxxxxyz[k];

                g_z_0_xxzzz_xxxxzz[k] = -g_z_0_xzzz_xxxxzz[k] * cd_x[k] + g_z_0_xzzz_xxxxxzz[k];

                g_z_0_xxzzz_xxxyyy[k] = -g_z_0_xzzz_xxxyyy[k] * cd_x[k] + g_z_0_xzzz_xxxxyyy[k];

                g_z_0_xxzzz_xxxyyz[k] = -g_z_0_xzzz_xxxyyz[k] * cd_x[k] + g_z_0_xzzz_xxxxyyz[k];

                g_z_0_xxzzz_xxxyzz[k] = -g_z_0_xzzz_xxxyzz[k] * cd_x[k] + g_z_0_xzzz_xxxxyzz[k];

                g_z_0_xxzzz_xxxzzz[k] = -g_z_0_xzzz_xxxzzz[k] * cd_x[k] + g_z_0_xzzz_xxxxzzz[k];

                g_z_0_xxzzz_xxyyyy[k] = -g_z_0_xzzz_xxyyyy[k] * cd_x[k] + g_z_0_xzzz_xxxyyyy[k];

                g_z_0_xxzzz_xxyyyz[k] = -g_z_0_xzzz_xxyyyz[k] * cd_x[k] + g_z_0_xzzz_xxxyyyz[k];

                g_z_0_xxzzz_xxyyzz[k] = -g_z_0_xzzz_xxyyzz[k] * cd_x[k] + g_z_0_xzzz_xxxyyzz[k];

                g_z_0_xxzzz_xxyzzz[k] = -g_z_0_xzzz_xxyzzz[k] * cd_x[k] + g_z_0_xzzz_xxxyzzz[k];

                g_z_0_xxzzz_xxzzzz[k] = -g_z_0_xzzz_xxzzzz[k] * cd_x[k] + g_z_0_xzzz_xxxzzzz[k];

                g_z_0_xxzzz_xyyyyy[k] = -g_z_0_xzzz_xyyyyy[k] * cd_x[k] + g_z_0_xzzz_xxyyyyy[k];

                g_z_0_xxzzz_xyyyyz[k] = -g_z_0_xzzz_xyyyyz[k] * cd_x[k] + g_z_0_xzzz_xxyyyyz[k];

                g_z_0_xxzzz_xyyyzz[k] = -g_z_0_xzzz_xyyyzz[k] * cd_x[k] + g_z_0_xzzz_xxyyyzz[k];

                g_z_0_xxzzz_xyyzzz[k] = -g_z_0_xzzz_xyyzzz[k] * cd_x[k] + g_z_0_xzzz_xxyyzzz[k];

                g_z_0_xxzzz_xyzzzz[k] = -g_z_0_xzzz_xyzzzz[k] * cd_x[k] + g_z_0_xzzz_xxyzzzz[k];

                g_z_0_xxzzz_xzzzzz[k] = -g_z_0_xzzz_xzzzzz[k] * cd_x[k] + g_z_0_xzzz_xxzzzzz[k];

                g_z_0_xxzzz_yyyyyy[k] = -g_z_0_xzzz_yyyyyy[k] * cd_x[k] + g_z_0_xzzz_xyyyyyy[k];

                g_z_0_xxzzz_yyyyyz[k] = -g_z_0_xzzz_yyyyyz[k] * cd_x[k] + g_z_0_xzzz_xyyyyyz[k];

                g_z_0_xxzzz_yyyyzz[k] = -g_z_0_xzzz_yyyyzz[k] * cd_x[k] + g_z_0_xzzz_xyyyyzz[k];

                g_z_0_xxzzz_yyyzzz[k] = -g_z_0_xzzz_yyyzzz[k] * cd_x[k] + g_z_0_xzzz_xyyyzzz[k];

                g_z_0_xxzzz_yyzzzz[k] = -g_z_0_xzzz_yyzzzz[k] * cd_x[k] + g_z_0_xzzz_xyyzzzz[k];

                g_z_0_xxzzz_yzzzzz[k] = -g_z_0_xzzz_yzzzzz[k] * cd_x[k] + g_z_0_xzzz_xyzzzzz[k];

                g_z_0_xxzzz_zzzzzz[k] = -g_z_0_xzzz_zzzzzz[k] * cd_x[k] + g_z_0_xzzz_xzzzzzz[k];
            }

            /// Set up 280-308 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 280);

            auto g_z_0_xyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 281);

            auto g_z_0_xyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 282);

            auto g_z_0_xyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 283);

            auto g_z_0_xyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 284);

            auto g_z_0_xyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 285);

            auto g_z_0_xyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 286);

            auto g_z_0_xyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 287);

            auto g_z_0_xyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 288);

            auto g_z_0_xyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 289);

            auto g_z_0_xyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 290);

            auto g_z_0_xyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 291);

            auto g_z_0_xyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 292);

            auto g_z_0_xyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 293);

            auto g_z_0_xyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 294);

            auto g_z_0_xyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 295);

            auto g_z_0_xyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 296);

            auto g_z_0_xyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 297);

            auto g_z_0_xyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 298);

            auto g_z_0_xyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 299);

            auto g_z_0_xyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 300);

            auto g_z_0_xyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 301);

            auto g_z_0_xyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 302);

            auto g_z_0_xyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 303);

            auto g_z_0_xyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 304);

            auto g_z_0_xyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 305);

            auto g_z_0_xyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 306);

            auto g_z_0_xyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 307);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyy_xxxxxx, g_z_0_xyyyy_xxxxxy, g_z_0_xyyyy_xxxxxz, g_z_0_xyyyy_xxxxyy, g_z_0_xyyyy_xxxxyz, g_z_0_xyyyy_xxxxzz, g_z_0_xyyyy_xxxyyy, g_z_0_xyyyy_xxxyyz, g_z_0_xyyyy_xxxyzz, g_z_0_xyyyy_xxxzzz, g_z_0_xyyyy_xxyyyy, g_z_0_xyyyy_xxyyyz, g_z_0_xyyyy_xxyyzz, g_z_0_xyyyy_xxyzzz, g_z_0_xyyyy_xxzzzz, g_z_0_xyyyy_xyyyyy, g_z_0_xyyyy_xyyyyz, g_z_0_xyyyy_xyyyzz, g_z_0_xyyyy_xyyzzz, g_z_0_xyyyy_xyzzzz, g_z_0_xyyyy_xzzzzz, g_z_0_xyyyy_yyyyyy, g_z_0_xyyyy_yyyyyz, g_z_0_xyyyy_yyyyzz, g_z_0_xyyyy_yyyzzz, g_z_0_xyyyy_yyzzzz, g_z_0_xyyyy_yzzzzz, g_z_0_xyyyy_zzzzzz, g_z_0_yyyy_xxxxxx, g_z_0_yyyy_xxxxxxx, g_z_0_yyyy_xxxxxxy, g_z_0_yyyy_xxxxxxz, g_z_0_yyyy_xxxxxy, g_z_0_yyyy_xxxxxyy, g_z_0_yyyy_xxxxxyz, g_z_0_yyyy_xxxxxz, g_z_0_yyyy_xxxxxzz, g_z_0_yyyy_xxxxyy, g_z_0_yyyy_xxxxyyy, g_z_0_yyyy_xxxxyyz, g_z_0_yyyy_xxxxyz, g_z_0_yyyy_xxxxyzz, g_z_0_yyyy_xxxxzz, g_z_0_yyyy_xxxxzzz, g_z_0_yyyy_xxxyyy, g_z_0_yyyy_xxxyyyy, g_z_0_yyyy_xxxyyyz, g_z_0_yyyy_xxxyyz, g_z_0_yyyy_xxxyyzz, g_z_0_yyyy_xxxyzz, g_z_0_yyyy_xxxyzzz, g_z_0_yyyy_xxxzzz, g_z_0_yyyy_xxxzzzz, g_z_0_yyyy_xxyyyy, g_z_0_yyyy_xxyyyyy, g_z_0_yyyy_xxyyyyz, g_z_0_yyyy_xxyyyz, g_z_0_yyyy_xxyyyzz, g_z_0_yyyy_xxyyzz, g_z_0_yyyy_xxyyzzz, g_z_0_yyyy_xxyzzz, g_z_0_yyyy_xxyzzzz, g_z_0_yyyy_xxzzzz, g_z_0_yyyy_xxzzzzz, g_z_0_yyyy_xyyyyy, g_z_0_yyyy_xyyyyyy, g_z_0_yyyy_xyyyyyz, g_z_0_yyyy_xyyyyz, g_z_0_yyyy_xyyyyzz, g_z_0_yyyy_xyyyzz, g_z_0_yyyy_xyyyzzz, g_z_0_yyyy_xyyzzz, g_z_0_yyyy_xyyzzzz, g_z_0_yyyy_xyzzzz, g_z_0_yyyy_xyzzzzz, g_z_0_yyyy_xzzzzz, g_z_0_yyyy_xzzzzzz, g_z_0_yyyy_yyyyyy, g_z_0_yyyy_yyyyyz, g_z_0_yyyy_yyyyzz, g_z_0_yyyy_yyyzzz, g_z_0_yyyy_yyzzzz, g_z_0_yyyy_yzzzzz, g_z_0_yyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyy_xxxxxx[k] = -g_z_0_yyyy_xxxxxx[k] * cd_x[k] + g_z_0_yyyy_xxxxxxx[k];

                g_z_0_xyyyy_xxxxxy[k] = -g_z_0_yyyy_xxxxxy[k] * cd_x[k] + g_z_0_yyyy_xxxxxxy[k];

                g_z_0_xyyyy_xxxxxz[k] = -g_z_0_yyyy_xxxxxz[k] * cd_x[k] + g_z_0_yyyy_xxxxxxz[k];

                g_z_0_xyyyy_xxxxyy[k] = -g_z_0_yyyy_xxxxyy[k] * cd_x[k] + g_z_0_yyyy_xxxxxyy[k];

                g_z_0_xyyyy_xxxxyz[k] = -g_z_0_yyyy_xxxxyz[k] * cd_x[k] + g_z_0_yyyy_xxxxxyz[k];

                g_z_0_xyyyy_xxxxzz[k] = -g_z_0_yyyy_xxxxzz[k] * cd_x[k] + g_z_0_yyyy_xxxxxzz[k];

                g_z_0_xyyyy_xxxyyy[k] = -g_z_0_yyyy_xxxyyy[k] * cd_x[k] + g_z_0_yyyy_xxxxyyy[k];

                g_z_0_xyyyy_xxxyyz[k] = -g_z_0_yyyy_xxxyyz[k] * cd_x[k] + g_z_0_yyyy_xxxxyyz[k];

                g_z_0_xyyyy_xxxyzz[k] = -g_z_0_yyyy_xxxyzz[k] * cd_x[k] + g_z_0_yyyy_xxxxyzz[k];

                g_z_0_xyyyy_xxxzzz[k] = -g_z_0_yyyy_xxxzzz[k] * cd_x[k] + g_z_0_yyyy_xxxxzzz[k];

                g_z_0_xyyyy_xxyyyy[k] = -g_z_0_yyyy_xxyyyy[k] * cd_x[k] + g_z_0_yyyy_xxxyyyy[k];

                g_z_0_xyyyy_xxyyyz[k] = -g_z_0_yyyy_xxyyyz[k] * cd_x[k] + g_z_0_yyyy_xxxyyyz[k];

                g_z_0_xyyyy_xxyyzz[k] = -g_z_0_yyyy_xxyyzz[k] * cd_x[k] + g_z_0_yyyy_xxxyyzz[k];

                g_z_0_xyyyy_xxyzzz[k] = -g_z_0_yyyy_xxyzzz[k] * cd_x[k] + g_z_0_yyyy_xxxyzzz[k];

                g_z_0_xyyyy_xxzzzz[k] = -g_z_0_yyyy_xxzzzz[k] * cd_x[k] + g_z_0_yyyy_xxxzzzz[k];

                g_z_0_xyyyy_xyyyyy[k] = -g_z_0_yyyy_xyyyyy[k] * cd_x[k] + g_z_0_yyyy_xxyyyyy[k];

                g_z_0_xyyyy_xyyyyz[k] = -g_z_0_yyyy_xyyyyz[k] * cd_x[k] + g_z_0_yyyy_xxyyyyz[k];

                g_z_0_xyyyy_xyyyzz[k] = -g_z_0_yyyy_xyyyzz[k] * cd_x[k] + g_z_0_yyyy_xxyyyzz[k];

                g_z_0_xyyyy_xyyzzz[k] = -g_z_0_yyyy_xyyzzz[k] * cd_x[k] + g_z_0_yyyy_xxyyzzz[k];

                g_z_0_xyyyy_xyzzzz[k] = -g_z_0_yyyy_xyzzzz[k] * cd_x[k] + g_z_0_yyyy_xxyzzzz[k];

                g_z_0_xyyyy_xzzzzz[k] = -g_z_0_yyyy_xzzzzz[k] * cd_x[k] + g_z_0_yyyy_xxzzzzz[k];

                g_z_0_xyyyy_yyyyyy[k] = -g_z_0_yyyy_yyyyyy[k] * cd_x[k] + g_z_0_yyyy_xyyyyyy[k];

                g_z_0_xyyyy_yyyyyz[k] = -g_z_0_yyyy_yyyyyz[k] * cd_x[k] + g_z_0_yyyy_xyyyyyz[k];

                g_z_0_xyyyy_yyyyzz[k] = -g_z_0_yyyy_yyyyzz[k] * cd_x[k] + g_z_0_yyyy_xyyyyzz[k];

                g_z_0_xyyyy_yyyzzz[k] = -g_z_0_yyyy_yyyzzz[k] * cd_x[k] + g_z_0_yyyy_xyyyzzz[k];

                g_z_0_xyyyy_yyzzzz[k] = -g_z_0_yyyy_yyzzzz[k] * cd_x[k] + g_z_0_yyyy_xyyzzzz[k];

                g_z_0_xyyyy_yzzzzz[k] = -g_z_0_yyyy_yzzzzz[k] * cd_x[k] + g_z_0_yyyy_xyzzzzz[k];

                g_z_0_xyyyy_zzzzzz[k] = -g_z_0_yyyy_zzzzzz[k] * cd_x[k] + g_z_0_yyyy_xzzzzzz[k];
            }

            /// Set up 308-336 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 308);

            auto g_z_0_xyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 309);

            auto g_z_0_xyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 310);

            auto g_z_0_xyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 311);

            auto g_z_0_xyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 312);

            auto g_z_0_xyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 313);

            auto g_z_0_xyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 314);

            auto g_z_0_xyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 315);

            auto g_z_0_xyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 316);

            auto g_z_0_xyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 317);

            auto g_z_0_xyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 318);

            auto g_z_0_xyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 319);

            auto g_z_0_xyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 320);

            auto g_z_0_xyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 321);

            auto g_z_0_xyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 322);

            auto g_z_0_xyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 323);

            auto g_z_0_xyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 324);

            auto g_z_0_xyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 325);

            auto g_z_0_xyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 326);

            auto g_z_0_xyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 327);

            auto g_z_0_xyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 328);

            auto g_z_0_xyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 329);

            auto g_z_0_xyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 330);

            auto g_z_0_xyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 331);

            auto g_z_0_xyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 332);

            auto g_z_0_xyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 333);

            auto g_z_0_xyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 334);

            auto g_z_0_xyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 335);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyz_xxxxxx, g_z_0_xyyyz_xxxxxy, g_z_0_xyyyz_xxxxxz, g_z_0_xyyyz_xxxxyy, g_z_0_xyyyz_xxxxyz, g_z_0_xyyyz_xxxxzz, g_z_0_xyyyz_xxxyyy, g_z_0_xyyyz_xxxyyz, g_z_0_xyyyz_xxxyzz, g_z_0_xyyyz_xxxzzz, g_z_0_xyyyz_xxyyyy, g_z_0_xyyyz_xxyyyz, g_z_0_xyyyz_xxyyzz, g_z_0_xyyyz_xxyzzz, g_z_0_xyyyz_xxzzzz, g_z_0_xyyyz_xyyyyy, g_z_0_xyyyz_xyyyyz, g_z_0_xyyyz_xyyyzz, g_z_0_xyyyz_xyyzzz, g_z_0_xyyyz_xyzzzz, g_z_0_xyyyz_xzzzzz, g_z_0_xyyyz_yyyyyy, g_z_0_xyyyz_yyyyyz, g_z_0_xyyyz_yyyyzz, g_z_0_xyyyz_yyyzzz, g_z_0_xyyyz_yyzzzz, g_z_0_xyyyz_yzzzzz, g_z_0_xyyyz_zzzzzz, g_z_0_yyyz_xxxxxx, g_z_0_yyyz_xxxxxxx, g_z_0_yyyz_xxxxxxy, g_z_0_yyyz_xxxxxxz, g_z_0_yyyz_xxxxxy, g_z_0_yyyz_xxxxxyy, g_z_0_yyyz_xxxxxyz, g_z_0_yyyz_xxxxxz, g_z_0_yyyz_xxxxxzz, g_z_0_yyyz_xxxxyy, g_z_0_yyyz_xxxxyyy, g_z_0_yyyz_xxxxyyz, g_z_0_yyyz_xxxxyz, g_z_0_yyyz_xxxxyzz, g_z_0_yyyz_xxxxzz, g_z_0_yyyz_xxxxzzz, g_z_0_yyyz_xxxyyy, g_z_0_yyyz_xxxyyyy, g_z_0_yyyz_xxxyyyz, g_z_0_yyyz_xxxyyz, g_z_0_yyyz_xxxyyzz, g_z_0_yyyz_xxxyzz, g_z_0_yyyz_xxxyzzz, g_z_0_yyyz_xxxzzz, g_z_0_yyyz_xxxzzzz, g_z_0_yyyz_xxyyyy, g_z_0_yyyz_xxyyyyy, g_z_0_yyyz_xxyyyyz, g_z_0_yyyz_xxyyyz, g_z_0_yyyz_xxyyyzz, g_z_0_yyyz_xxyyzz, g_z_0_yyyz_xxyyzzz, g_z_0_yyyz_xxyzzz, g_z_0_yyyz_xxyzzzz, g_z_0_yyyz_xxzzzz, g_z_0_yyyz_xxzzzzz, g_z_0_yyyz_xyyyyy, g_z_0_yyyz_xyyyyyy, g_z_0_yyyz_xyyyyyz, g_z_0_yyyz_xyyyyz, g_z_0_yyyz_xyyyyzz, g_z_0_yyyz_xyyyzz, g_z_0_yyyz_xyyyzzz, g_z_0_yyyz_xyyzzz, g_z_0_yyyz_xyyzzzz, g_z_0_yyyz_xyzzzz, g_z_0_yyyz_xyzzzzz, g_z_0_yyyz_xzzzzz, g_z_0_yyyz_xzzzzzz, g_z_0_yyyz_yyyyyy, g_z_0_yyyz_yyyyyz, g_z_0_yyyz_yyyyzz, g_z_0_yyyz_yyyzzz, g_z_0_yyyz_yyzzzz, g_z_0_yyyz_yzzzzz, g_z_0_yyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyz_xxxxxx[k] = -g_z_0_yyyz_xxxxxx[k] * cd_x[k] + g_z_0_yyyz_xxxxxxx[k];

                g_z_0_xyyyz_xxxxxy[k] = -g_z_0_yyyz_xxxxxy[k] * cd_x[k] + g_z_0_yyyz_xxxxxxy[k];

                g_z_0_xyyyz_xxxxxz[k] = -g_z_0_yyyz_xxxxxz[k] * cd_x[k] + g_z_0_yyyz_xxxxxxz[k];

                g_z_0_xyyyz_xxxxyy[k] = -g_z_0_yyyz_xxxxyy[k] * cd_x[k] + g_z_0_yyyz_xxxxxyy[k];

                g_z_0_xyyyz_xxxxyz[k] = -g_z_0_yyyz_xxxxyz[k] * cd_x[k] + g_z_0_yyyz_xxxxxyz[k];

                g_z_0_xyyyz_xxxxzz[k] = -g_z_0_yyyz_xxxxzz[k] * cd_x[k] + g_z_0_yyyz_xxxxxzz[k];

                g_z_0_xyyyz_xxxyyy[k] = -g_z_0_yyyz_xxxyyy[k] * cd_x[k] + g_z_0_yyyz_xxxxyyy[k];

                g_z_0_xyyyz_xxxyyz[k] = -g_z_0_yyyz_xxxyyz[k] * cd_x[k] + g_z_0_yyyz_xxxxyyz[k];

                g_z_0_xyyyz_xxxyzz[k] = -g_z_0_yyyz_xxxyzz[k] * cd_x[k] + g_z_0_yyyz_xxxxyzz[k];

                g_z_0_xyyyz_xxxzzz[k] = -g_z_0_yyyz_xxxzzz[k] * cd_x[k] + g_z_0_yyyz_xxxxzzz[k];

                g_z_0_xyyyz_xxyyyy[k] = -g_z_0_yyyz_xxyyyy[k] * cd_x[k] + g_z_0_yyyz_xxxyyyy[k];

                g_z_0_xyyyz_xxyyyz[k] = -g_z_0_yyyz_xxyyyz[k] * cd_x[k] + g_z_0_yyyz_xxxyyyz[k];

                g_z_0_xyyyz_xxyyzz[k] = -g_z_0_yyyz_xxyyzz[k] * cd_x[k] + g_z_0_yyyz_xxxyyzz[k];

                g_z_0_xyyyz_xxyzzz[k] = -g_z_0_yyyz_xxyzzz[k] * cd_x[k] + g_z_0_yyyz_xxxyzzz[k];

                g_z_0_xyyyz_xxzzzz[k] = -g_z_0_yyyz_xxzzzz[k] * cd_x[k] + g_z_0_yyyz_xxxzzzz[k];

                g_z_0_xyyyz_xyyyyy[k] = -g_z_0_yyyz_xyyyyy[k] * cd_x[k] + g_z_0_yyyz_xxyyyyy[k];

                g_z_0_xyyyz_xyyyyz[k] = -g_z_0_yyyz_xyyyyz[k] * cd_x[k] + g_z_0_yyyz_xxyyyyz[k];

                g_z_0_xyyyz_xyyyzz[k] = -g_z_0_yyyz_xyyyzz[k] * cd_x[k] + g_z_0_yyyz_xxyyyzz[k];

                g_z_0_xyyyz_xyyzzz[k] = -g_z_0_yyyz_xyyzzz[k] * cd_x[k] + g_z_0_yyyz_xxyyzzz[k];

                g_z_0_xyyyz_xyzzzz[k] = -g_z_0_yyyz_xyzzzz[k] * cd_x[k] + g_z_0_yyyz_xxyzzzz[k];

                g_z_0_xyyyz_xzzzzz[k] = -g_z_0_yyyz_xzzzzz[k] * cd_x[k] + g_z_0_yyyz_xxzzzzz[k];

                g_z_0_xyyyz_yyyyyy[k] = -g_z_0_yyyz_yyyyyy[k] * cd_x[k] + g_z_0_yyyz_xyyyyyy[k];

                g_z_0_xyyyz_yyyyyz[k] = -g_z_0_yyyz_yyyyyz[k] * cd_x[k] + g_z_0_yyyz_xyyyyyz[k];

                g_z_0_xyyyz_yyyyzz[k] = -g_z_0_yyyz_yyyyzz[k] * cd_x[k] + g_z_0_yyyz_xyyyyzz[k];

                g_z_0_xyyyz_yyyzzz[k] = -g_z_0_yyyz_yyyzzz[k] * cd_x[k] + g_z_0_yyyz_xyyyzzz[k];

                g_z_0_xyyyz_yyzzzz[k] = -g_z_0_yyyz_yyzzzz[k] * cd_x[k] + g_z_0_yyyz_xyyzzzz[k];

                g_z_0_xyyyz_yzzzzz[k] = -g_z_0_yyyz_yzzzzz[k] * cd_x[k] + g_z_0_yyyz_xyzzzzz[k];

                g_z_0_xyyyz_zzzzzz[k] = -g_z_0_yyyz_zzzzzz[k] * cd_x[k] + g_z_0_yyyz_xzzzzzz[k];
            }

            /// Set up 336-364 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 336);

            auto g_z_0_xyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 337);

            auto g_z_0_xyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 338);

            auto g_z_0_xyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 339);

            auto g_z_0_xyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 340);

            auto g_z_0_xyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 341);

            auto g_z_0_xyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 342);

            auto g_z_0_xyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 343);

            auto g_z_0_xyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 344);

            auto g_z_0_xyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 345);

            auto g_z_0_xyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 346);

            auto g_z_0_xyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 347);

            auto g_z_0_xyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 348);

            auto g_z_0_xyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 349);

            auto g_z_0_xyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 350);

            auto g_z_0_xyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 351);

            auto g_z_0_xyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 352);

            auto g_z_0_xyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 353);

            auto g_z_0_xyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 354);

            auto g_z_0_xyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 355);

            auto g_z_0_xyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 356);

            auto g_z_0_xyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 357);

            auto g_z_0_xyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 358);

            auto g_z_0_xyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 359);

            auto g_z_0_xyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 360);

            auto g_z_0_xyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 361);

            auto g_z_0_xyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 362);

            auto g_z_0_xyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 363);

            #pragma omp simd aligned(cd_x, g_z_0_xyyzz_xxxxxx, g_z_0_xyyzz_xxxxxy, g_z_0_xyyzz_xxxxxz, g_z_0_xyyzz_xxxxyy, g_z_0_xyyzz_xxxxyz, g_z_0_xyyzz_xxxxzz, g_z_0_xyyzz_xxxyyy, g_z_0_xyyzz_xxxyyz, g_z_0_xyyzz_xxxyzz, g_z_0_xyyzz_xxxzzz, g_z_0_xyyzz_xxyyyy, g_z_0_xyyzz_xxyyyz, g_z_0_xyyzz_xxyyzz, g_z_0_xyyzz_xxyzzz, g_z_0_xyyzz_xxzzzz, g_z_0_xyyzz_xyyyyy, g_z_0_xyyzz_xyyyyz, g_z_0_xyyzz_xyyyzz, g_z_0_xyyzz_xyyzzz, g_z_0_xyyzz_xyzzzz, g_z_0_xyyzz_xzzzzz, g_z_0_xyyzz_yyyyyy, g_z_0_xyyzz_yyyyyz, g_z_0_xyyzz_yyyyzz, g_z_0_xyyzz_yyyzzz, g_z_0_xyyzz_yyzzzz, g_z_0_xyyzz_yzzzzz, g_z_0_xyyzz_zzzzzz, g_z_0_yyzz_xxxxxx, g_z_0_yyzz_xxxxxxx, g_z_0_yyzz_xxxxxxy, g_z_0_yyzz_xxxxxxz, g_z_0_yyzz_xxxxxy, g_z_0_yyzz_xxxxxyy, g_z_0_yyzz_xxxxxyz, g_z_0_yyzz_xxxxxz, g_z_0_yyzz_xxxxxzz, g_z_0_yyzz_xxxxyy, g_z_0_yyzz_xxxxyyy, g_z_0_yyzz_xxxxyyz, g_z_0_yyzz_xxxxyz, g_z_0_yyzz_xxxxyzz, g_z_0_yyzz_xxxxzz, g_z_0_yyzz_xxxxzzz, g_z_0_yyzz_xxxyyy, g_z_0_yyzz_xxxyyyy, g_z_0_yyzz_xxxyyyz, g_z_0_yyzz_xxxyyz, g_z_0_yyzz_xxxyyzz, g_z_0_yyzz_xxxyzz, g_z_0_yyzz_xxxyzzz, g_z_0_yyzz_xxxzzz, g_z_0_yyzz_xxxzzzz, g_z_0_yyzz_xxyyyy, g_z_0_yyzz_xxyyyyy, g_z_0_yyzz_xxyyyyz, g_z_0_yyzz_xxyyyz, g_z_0_yyzz_xxyyyzz, g_z_0_yyzz_xxyyzz, g_z_0_yyzz_xxyyzzz, g_z_0_yyzz_xxyzzz, g_z_0_yyzz_xxyzzzz, g_z_0_yyzz_xxzzzz, g_z_0_yyzz_xxzzzzz, g_z_0_yyzz_xyyyyy, g_z_0_yyzz_xyyyyyy, g_z_0_yyzz_xyyyyyz, g_z_0_yyzz_xyyyyz, g_z_0_yyzz_xyyyyzz, g_z_0_yyzz_xyyyzz, g_z_0_yyzz_xyyyzzz, g_z_0_yyzz_xyyzzz, g_z_0_yyzz_xyyzzzz, g_z_0_yyzz_xyzzzz, g_z_0_yyzz_xyzzzzz, g_z_0_yyzz_xzzzzz, g_z_0_yyzz_xzzzzzz, g_z_0_yyzz_yyyyyy, g_z_0_yyzz_yyyyyz, g_z_0_yyzz_yyyyzz, g_z_0_yyzz_yyyzzz, g_z_0_yyzz_yyzzzz, g_z_0_yyzz_yzzzzz, g_z_0_yyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzz_xxxxxx[k] = -g_z_0_yyzz_xxxxxx[k] * cd_x[k] + g_z_0_yyzz_xxxxxxx[k];

                g_z_0_xyyzz_xxxxxy[k] = -g_z_0_yyzz_xxxxxy[k] * cd_x[k] + g_z_0_yyzz_xxxxxxy[k];

                g_z_0_xyyzz_xxxxxz[k] = -g_z_0_yyzz_xxxxxz[k] * cd_x[k] + g_z_0_yyzz_xxxxxxz[k];

                g_z_0_xyyzz_xxxxyy[k] = -g_z_0_yyzz_xxxxyy[k] * cd_x[k] + g_z_0_yyzz_xxxxxyy[k];

                g_z_0_xyyzz_xxxxyz[k] = -g_z_0_yyzz_xxxxyz[k] * cd_x[k] + g_z_0_yyzz_xxxxxyz[k];

                g_z_0_xyyzz_xxxxzz[k] = -g_z_0_yyzz_xxxxzz[k] * cd_x[k] + g_z_0_yyzz_xxxxxzz[k];

                g_z_0_xyyzz_xxxyyy[k] = -g_z_0_yyzz_xxxyyy[k] * cd_x[k] + g_z_0_yyzz_xxxxyyy[k];

                g_z_0_xyyzz_xxxyyz[k] = -g_z_0_yyzz_xxxyyz[k] * cd_x[k] + g_z_0_yyzz_xxxxyyz[k];

                g_z_0_xyyzz_xxxyzz[k] = -g_z_0_yyzz_xxxyzz[k] * cd_x[k] + g_z_0_yyzz_xxxxyzz[k];

                g_z_0_xyyzz_xxxzzz[k] = -g_z_0_yyzz_xxxzzz[k] * cd_x[k] + g_z_0_yyzz_xxxxzzz[k];

                g_z_0_xyyzz_xxyyyy[k] = -g_z_0_yyzz_xxyyyy[k] * cd_x[k] + g_z_0_yyzz_xxxyyyy[k];

                g_z_0_xyyzz_xxyyyz[k] = -g_z_0_yyzz_xxyyyz[k] * cd_x[k] + g_z_0_yyzz_xxxyyyz[k];

                g_z_0_xyyzz_xxyyzz[k] = -g_z_0_yyzz_xxyyzz[k] * cd_x[k] + g_z_0_yyzz_xxxyyzz[k];

                g_z_0_xyyzz_xxyzzz[k] = -g_z_0_yyzz_xxyzzz[k] * cd_x[k] + g_z_0_yyzz_xxxyzzz[k];

                g_z_0_xyyzz_xxzzzz[k] = -g_z_0_yyzz_xxzzzz[k] * cd_x[k] + g_z_0_yyzz_xxxzzzz[k];

                g_z_0_xyyzz_xyyyyy[k] = -g_z_0_yyzz_xyyyyy[k] * cd_x[k] + g_z_0_yyzz_xxyyyyy[k];

                g_z_0_xyyzz_xyyyyz[k] = -g_z_0_yyzz_xyyyyz[k] * cd_x[k] + g_z_0_yyzz_xxyyyyz[k];

                g_z_0_xyyzz_xyyyzz[k] = -g_z_0_yyzz_xyyyzz[k] * cd_x[k] + g_z_0_yyzz_xxyyyzz[k];

                g_z_0_xyyzz_xyyzzz[k] = -g_z_0_yyzz_xyyzzz[k] * cd_x[k] + g_z_0_yyzz_xxyyzzz[k];

                g_z_0_xyyzz_xyzzzz[k] = -g_z_0_yyzz_xyzzzz[k] * cd_x[k] + g_z_0_yyzz_xxyzzzz[k];

                g_z_0_xyyzz_xzzzzz[k] = -g_z_0_yyzz_xzzzzz[k] * cd_x[k] + g_z_0_yyzz_xxzzzzz[k];

                g_z_0_xyyzz_yyyyyy[k] = -g_z_0_yyzz_yyyyyy[k] * cd_x[k] + g_z_0_yyzz_xyyyyyy[k];

                g_z_0_xyyzz_yyyyyz[k] = -g_z_0_yyzz_yyyyyz[k] * cd_x[k] + g_z_0_yyzz_xyyyyyz[k];

                g_z_0_xyyzz_yyyyzz[k] = -g_z_0_yyzz_yyyyzz[k] * cd_x[k] + g_z_0_yyzz_xyyyyzz[k];

                g_z_0_xyyzz_yyyzzz[k] = -g_z_0_yyzz_yyyzzz[k] * cd_x[k] + g_z_0_yyzz_xyyyzzz[k];

                g_z_0_xyyzz_yyzzzz[k] = -g_z_0_yyzz_yyzzzz[k] * cd_x[k] + g_z_0_yyzz_xyyzzzz[k];

                g_z_0_xyyzz_yzzzzz[k] = -g_z_0_yyzz_yzzzzz[k] * cd_x[k] + g_z_0_yyzz_xyzzzzz[k];

                g_z_0_xyyzz_zzzzzz[k] = -g_z_0_yyzz_zzzzzz[k] * cd_x[k] + g_z_0_yyzz_xzzzzzz[k];
            }

            /// Set up 364-392 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 364);

            auto g_z_0_xyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 365);

            auto g_z_0_xyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 366);

            auto g_z_0_xyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 367);

            auto g_z_0_xyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 368);

            auto g_z_0_xyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 369);

            auto g_z_0_xyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 370);

            auto g_z_0_xyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 371);

            auto g_z_0_xyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 372);

            auto g_z_0_xyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 373);

            auto g_z_0_xyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 374);

            auto g_z_0_xyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 375);

            auto g_z_0_xyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 376);

            auto g_z_0_xyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 377);

            auto g_z_0_xyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 378);

            auto g_z_0_xyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 379);

            auto g_z_0_xyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 380);

            auto g_z_0_xyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 381);

            auto g_z_0_xyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 382);

            auto g_z_0_xyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 383);

            auto g_z_0_xyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 384);

            auto g_z_0_xyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 385);

            auto g_z_0_xyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 386);

            auto g_z_0_xyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 387);

            auto g_z_0_xyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 388);

            auto g_z_0_xyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 389);

            auto g_z_0_xyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 390);

            auto g_z_0_xyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 391);

            #pragma omp simd aligned(cd_x, g_z_0_xyzzz_xxxxxx, g_z_0_xyzzz_xxxxxy, g_z_0_xyzzz_xxxxxz, g_z_0_xyzzz_xxxxyy, g_z_0_xyzzz_xxxxyz, g_z_0_xyzzz_xxxxzz, g_z_0_xyzzz_xxxyyy, g_z_0_xyzzz_xxxyyz, g_z_0_xyzzz_xxxyzz, g_z_0_xyzzz_xxxzzz, g_z_0_xyzzz_xxyyyy, g_z_0_xyzzz_xxyyyz, g_z_0_xyzzz_xxyyzz, g_z_0_xyzzz_xxyzzz, g_z_0_xyzzz_xxzzzz, g_z_0_xyzzz_xyyyyy, g_z_0_xyzzz_xyyyyz, g_z_0_xyzzz_xyyyzz, g_z_0_xyzzz_xyyzzz, g_z_0_xyzzz_xyzzzz, g_z_0_xyzzz_xzzzzz, g_z_0_xyzzz_yyyyyy, g_z_0_xyzzz_yyyyyz, g_z_0_xyzzz_yyyyzz, g_z_0_xyzzz_yyyzzz, g_z_0_xyzzz_yyzzzz, g_z_0_xyzzz_yzzzzz, g_z_0_xyzzz_zzzzzz, g_z_0_yzzz_xxxxxx, g_z_0_yzzz_xxxxxxx, g_z_0_yzzz_xxxxxxy, g_z_0_yzzz_xxxxxxz, g_z_0_yzzz_xxxxxy, g_z_0_yzzz_xxxxxyy, g_z_0_yzzz_xxxxxyz, g_z_0_yzzz_xxxxxz, g_z_0_yzzz_xxxxxzz, g_z_0_yzzz_xxxxyy, g_z_0_yzzz_xxxxyyy, g_z_0_yzzz_xxxxyyz, g_z_0_yzzz_xxxxyz, g_z_0_yzzz_xxxxyzz, g_z_0_yzzz_xxxxzz, g_z_0_yzzz_xxxxzzz, g_z_0_yzzz_xxxyyy, g_z_0_yzzz_xxxyyyy, g_z_0_yzzz_xxxyyyz, g_z_0_yzzz_xxxyyz, g_z_0_yzzz_xxxyyzz, g_z_0_yzzz_xxxyzz, g_z_0_yzzz_xxxyzzz, g_z_0_yzzz_xxxzzz, g_z_0_yzzz_xxxzzzz, g_z_0_yzzz_xxyyyy, g_z_0_yzzz_xxyyyyy, g_z_0_yzzz_xxyyyyz, g_z_0_yzzz_xxyyyz, g_z_0_yzzz_xxyyyzz, g_z_0_yzzz_xxyyzz, g_z_0_yzzz_xxyyzzz, g_z_0_yzzz_xxyzzz, g_z_0_yzzz_xxyzzzz, g_z_0_yzzz_xxzzzz, g_z_0_yzzz_xxzzzzz, g_z_0_yzzz_xyyyyy, g_z_0_yzzz_xyyyyyy, g_z_0_yzzz_xyyyyyz, g_z_0_yzzz_xyyyyz, g_z_0_yzzz_xyyyyzz, g_z_0_yzzz_xyyyzz, g_z_0_yzzz_xyyyzzz, g_z_0_yzzz_xyyzzz, g_z_0_yzzz_xyyzzzz, g_z_0_yzzz_xyzzzz, g_z_0_yzzz_xyzzzzz, g_z_0_yzzz_xzzzzz, g_z_0_yzzz_xzzzzzz, g_z_0_yzzz_yyyyyy, g_z_0_yzzz_yyyyyz, g_z_0_yzzz_yyyyzz, g_z_0_yzzz_yyyzzz, g_z_0_yzzz_yyzzzz, g_z_0_yzzz_yzzzzz, g_z_0_yzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzz_xxxxxx[k] = -g_z_0_yzzz_xxxxxx[k] * cd_x[k] + g_z_0_yzzz_xxxxxxx[k];

                g_z_0_xyzzz_xxxxxy[k] = -g_z_0_yzzz_xxxxxy[k] * cd_x[k] + g_z_0_yzzz_xxxxxxy[k];

                g_z_0_xyzzz_xxxxxz[k] = -g_z_0_yzzz_xxxxxz[k] * cd_x[k] + g_z_0_yzzz_xxxxxxz[k];

                g_z_0_xyzzz_xxxxyy[k] = -g_z_0_yzzz_xxxxyy[k] * cd_x[k] + g_z_0_yzzz_xxxxxyy[k];

                g_z_0_xyzzz_xxxxyz[k] = -g_z_0_yzzz_xxxxyz[k] * cd_x[k] + g_z_0_yzzz_xxxxxyz[k];

                g_z_0_xyzzz_xxxxzz[k] = -g_z_0_yzzz_xxxxzz[k] * cd_x[k] + g_z_0_yzzz_xxxxxzz[k];

                g_z_0_xyzzz_xxxyyy[k] = -g_z_0_yzzz_xxxyyy[k] * cd_x[k] + g_z_0_yzzz_xxxxyyy[k];

                g_z_0_xyzzz_xxxyyz[k] = -g_z_0_yzzz_xxxyyz[k] * cd_x[k] + g_z_0_yzzz_xxxxyyz[k];

                g_z_0_xyzzz_xxxyzz[k] = -g_z_0_yzzz_xxxyzz[k] * cd_x[k] + g_z_0_yzzz_xxxxyzz[k];

                g_z_0_xyzzz_xxxzzz[k] = -g_z_0_yzzz_xxxzzz[k] * cd_x[k] + g_z_0_yzzz_xxxxzzz[k];

                g_z_0_xyzzz_xxyyyy[k] = -g_z_0_yzzz_xxyyyy[k] * cd_x[k] + g_z_0_yzzz_xxxyyyy[k];

                g_z_0_xyzzz_xxyyyz[k] = -g_z_0_yzzz_xxyyyz[k] * cd_x[k] + g_z_0_yzzz_xxxyyyz[k];

                g_z_0_xyzzz_xxyyzz[k] = -g_z_0_yzzz_xxyyzz[k] * cd_x[k] + g_z_0_yzzz_xxxyyzz[k];

                g_z_0_xyzzz_xxyzzz[k] = -g_z_0_yzzz_xxyzzz[k] * cd_x[k] + g_z_0_yzzz_xxxyzzz[k];

                g_z_0_xyzzz_xxzzzz[k] = -g_z_0_yzzz_xxzzzz[k] * cd_x[k] + g_z_0_yzzz_xxxzzzz[k];

                g_z_0_xyzzz_xyyyyy[k] = -g_z_0_yzzz_xyyyyy[k] * cd_x[k] + g_z_0_yzzz_xxyyyyy[k];

                g_z_0_xyzzz_xyyyyz[k] = -g_z_0_yzzz_xyyyyz[k] * cd_x[k] + g_z_0_yzzz_xxyyyyz[k];

                g_z_0_xyzzz_xyyyzz[k] = -g_z_0_yzzz_xyyyzz[k] * cd_x[k] + g_z_0_yzzz_xxyyyzz[k];

                g_z_0_xyzzz_xyyzzz[k] = -g_z_0_yzzz_xyyzzz[k] * cd_x[k] + g_z_0_yzzz_xxyyzzz[k];

                g_z_0_xyzzz_xyzzzz[k] = -g_z_0_yzzz_xyzzzz[k] * cd_x[k] + g_z_0_yzzz_xxyzzzz[k];

                g_z_0_xyzzz_xzzzzz[k] = -g_z_0_yzzz_xzzzzz[k] * cd_x[k] + g_z_0_yzzz_xxzzzzz[k];

                g_z_0_xyzzz_yyyyyy[k] = -g_z_0_yzzz_yyyyyy[k] * cd_x[k] + g_z_0_yzzz_xyyyyyy[k];

                g_z_0_xyzzz_yyyyyz[k] = -g_z_0_yzzz_yyyyyz[k] * cd_x[k] + g_z_0_yzzz_xyyyyyz[k];

                g_z_0_xyzzz_yyyyzz[k] = -g_z_0_yzzz_yyyyzz[k] * cd_x[k] + g_z_0_yzzz_xyyyyzz[k];

                g_z_0_xyzzz_yyyzzz[k] = -g_z_0_yzzz_yyyzzz[k] * cd_x[k] + g_z_0_yzzz_xyyyzzz[k];

                g_z_0_xyzzz_yyzzzz[k] = -g_z_0_yzzz_yyzzzz[k] * cd_x[k] + g_z_0_yzzz_xyyzzzz[k];

                g_z_0_xyzzz_yzzzzz[k] = -g_z_0_yzzz_yzzzzz[k] * cd_x[k] + g_z_0_yzzz_xyzzzzz[k];

                g_z_0_xyzzz_zzzzzz[k] = -g_z_0_yzzz_zzzzzz[k] * cd_x[k] + g_z_0_yzzz_xzzzzzz[k];
            }

            /// Set up 392-420 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 392);

            auto g_z_0_xzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 393);

            auto g_z_0_xzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 394);

            auto g_z_0_xzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 395);

            auto g_z_0_xzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 396);

            auto g_z_0_xzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 397);

            auto g_z_0_xzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 398);

            auto g_z_0_xzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 399);

            auto g_z_0_xzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 400);

            auto g_z_0_xzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 401);

            auto g_z_0_xzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 402);

            auto g_z_0_xzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 403);

            auto g_z_0_xzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 404);

            auto g_z_0_xzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 405);

            auto g_z_0_xzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 406);

            auto g_z_0_xzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 407);

            auto g_z_0_xzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 408);

            auto g_z_0_xzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 409);

            auto g_z_0_xzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 410);

            auto g_z_0_xzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 411);

            auto g_z_0_xzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 412);

            auto g_z_0_xzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 413);

            auto g_z_0_xzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 414);

            auto g_z_0_xzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 415);

            auto g_z_0_xzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 416);

            auto g_z_0_xzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 417);

            auto g_z_0_xzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 418);

            auto g_z_0_xzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 419);

            #pragma omp simd aligned(cd_x, g_z_0_xzzzz_xxxxxx, g_z_0_xzzzz_xxxxxy, g_z_0_xzzzz_xxxxxz, g_z_0_xzzzz_xxxxyy, g_z_0_xzzzz_xxxxyz, g_z_0_xzzzz_xxxxzz, g_z_0_xzzzz_xxxyyy, g_z_0_xzzzz_xxxyyz, g_z_0_xzzzz_xxxyzz, g_z_0_xzzzz_xxxzzz, g_z_0_xzzzz_xxyyyy, g_z_0_xzzzz_xxyyyz, g_z_0_xzzzz_xxyyzz, g_z_0_xzzzz_xxyzzz, g_z_0_xzzzz_xxzzzz, g_z_0_xzzzz_xyyyyy, g_z_0_xzzzz_xyyyyz, g_z_0_xzzzz_xyyyzz, g_z_0_xzzzz_xyyzzz, g_z_0_xzzzz_xyzzzz, g_z_0_xzzzz_xzzzzz, g_z_0_xzzzz_yyyyyy, g_z_0_xzzzz_yyyyyz, g_z_0_xzzzz_yyyyzz, g_z_0_xzzzz_yyyzzz, g_z_0_xzzzz_yyzzzz, g_z_0_xzzzz_yzzzzz, g_z_0_xzzzz_zzzzzz, g_z_0_zzzz_xxxxxx, g_z_0_zzzz_xxxxxxx, g_z_0_zzzz_xxxxxxy, g_z_0_zzzz_xxxxxxz, g_z_0_zzzz_xxxxxy, g_z_0_zzzz_xxxxxyy, g_z_0_zzzz_xxxxxyz, g_z_0_zzzz_xxxxxz, g_z_0_zzzz_xxxxxzz, g_z_0_zzzz_xxxxyy, g_z_0_zzzz_xxxxyyy, g_z_0_zzzz_xxxxyyz, g_z_0_zzzz_xxxxyz, g_z_0_zzzz_xxxxyzz, g_z_0_zzzz_xxxxzz, g_z_0_zzzz_xxxxzzz, g_z_0_zzzz_xxxyyy, g_z_0_zzzz_xxxyyyy, g_z_0_zzzz_xxxyyyz, g_z_0_zzzz_xxxyyz, g_z_0_zzzz_xxxyyzz, g_z_0_zzzz_xxxyzz, g_z_0_zzzz_xxxyzzz, g_z_0_zzzz_xxxzzz, g_z_0_zzzz_xxxzzzz, g_z_0_zzzz_xxyyyy, g_z_0_zzzz_xxyyyyy, g_z_0_zzzz_xxyyyyz, g_z_0_zzzz_xxyyyz, g_z_0_zzzz_xxyyyzz, g_z_0_zzzz_xxyyzz, g_z_0_zzzz_xxyyzzz, g_z_0_zzzz_xxyzzz, g_z_0_zzzz_xxyzzzz, g_z_0_zzzz_xxzzzz, g_z_0_zzzz_xxzzzzz, g_z_0_zzzz_xyyyyy, g_z_0_zzzz_xyyyyyy, g_z_0_zzzz_xyyyyyz, g_z_0_zzzz_xyyyyz, g_z_0_zzzz_xyyyyzz, g_z_0_zzzz_xyyyzz, g_z_0_zzzz_xyyyzzz, g_z_0_zzzz_xyyzzz, g_z_0_zzzz_xyyzzzz, g_z_0_zzzz_xyzzzz, g_z_0_zzzz_xyzzzzz, g_z_0_zzzz_xzzzzz, g_z_0_zzzz_xzzzzzz, g_z_0_zzzz_yyyyyy, g_z_0_zzzz_yyyyyz, g_z_0_zzzz_yyyyzz, g_z_0_zzzz_yyyzzz, g_z_0_zzzz_yyzzzz, g_z_0_zzzz_yzzzzz, g_z_0_zzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzz_xxxxxx[k] = -g_z_0_zzzz_xxxxxx[k] * cd_x[k] + g_z_0_zzzz_xxxxxxx[k];

                g_z_0_xzzzz_xxxxxy[k] = -g_z_0_zzzz_xxxxxy[k] * cd_x[k] + g_z_0_zzzz_xxxxxxy[k];

                g_z_0_xzzzz_xxxxxz[k] = -g_z_0_zzzz_xxxxxz[k] * cd_x[k] + g_z_0_zzzz_xxxxxxz[k];

                g_z_0_xzzzz_xxxxyy[k] = -g_z_0_zzzz_xxxxyy[k] * cd_x[k] + g_z_0_zzzz_xxxxxyy[k];

                g_z_0_xzzzz_xxxxyz[k] = -g_z_0_zzzz_xxxxyz[k] * cd_x[k] + g_z_0_zzzz_xxxxxyz[k];

                g_z_0_xzzzz_xxxxzz[k] = -g_z_0_zzzz_xxxxzz[k] * cd_x[k] + g_z_0_zzzz_xxxxxzz[k];

                g_z_0_xzzzz_xxxyyy[k] = -g_z_0_zzzz_xxxyyy[k] * cd_x[k] + g_z_0_zzzz_xxxxyyy[k];

                g_z_0_xzzzz_xxxyyz[k] = -g_z_0_zzzz_xxxyyz[k] * cd_x[k] + g_z_0_zzzz_xxxxyyz[k];

                g_z_0_xzzzz_xxxyzz[k] = -g_z_0_zzzz_xxxyzz[k] * cd_x[k] + g_z_0_zzzz_xxxxyzz[k];

                g_z_0_xzzzz_xxxzzz[k] = -g_z_0_zzzz_xxxzzz[k] * cd_x[k] + g_z_0_zzzz_xxxxzzz[k];

                g_z_0_xzzzz_xxyyyy[k] = -g_z_0_zzzz_xxyyyy[k] * cd_x[k] + g_z_0_zzzz_xxxyyyy[k];

                g_z_0_xzzzz_xxyyyz[k] = -g_z_0_zzzz_xxyyyz[k] * cd_x[k] + g_z_0_zzzz_xxxyyyz[k];

                g_z_0_xzzzz_xxyyzz[k] = -g_z_0_zzzz_xxyyzz[k] * cd_x[k] + g_z_0_zzzz_xxxyyzz[k];

                g_z_0_xzzzz_xxyzzz[k] = -g_z_0_zzzz_xxyzzz[k] * cd_x[k] + g_z_0_zzzz_xxxyzzz[k];

                g_z_0_xzzzz_xxzzzz[k] = -g_z_0_zzzz_xxzzzz[k] * cd_x[k] + g_z_0_zzzz_xxxzzzz[k];

                g_z_0_xzzzz_xyyyyy[k] = -g_z_0_zzzz_xyyyyy[k] * cd_x[k] + g_z_0_zzzz_xxyyyyy[k];

                g_z_0_xzzzz_xyyyyz[k] = -g_z_0_zzzz_xyyyyz[k] * cd_x[k] + g_z_0_zzzz_xxyyyyz[k];

                g_z_0_xzzzz_xyyyzz[k] = -g_z_0_zzzz_xyyyzz[k] * cd_x[k] + g_z_0_zzzz_xxyyyzz[k];

                g_z_0_xzzzz_xyyzzz[k] = -g_z_0_zzzz_xyyzzz[k] * cd_x[k] + g_z_0_zzzz_xxyyzzz[k];

                g_z_0_xzzzz_xyzzzz[k] = -g_z_0_zzzz_xyzzzz[k] * cd_x[k] + g_z_0_zzzz_xxyzzzz[k];

                g_z_0_xzzzz_xzzzzz[k] = -g_z_0_zzzz_xzzzzz[k] * cd_x[k] + g_z_0_zzzz_xxzzzzz[k];

                g_z_0_xzzzz_yyyyyy[k] = -g_z_0_zzzz_yyyyyy[k] * cd_x[k] + g_z_0_zzzz_xyyyyyy[k];

                g_z_0_xzzzz_yyyyyz[k] = -g_z_0_zzzz_yyyyyz[k] * cd_x[k] + g_z_0_zzzz_xyyyyyz[k];

                g_z_0_xzzzz_yyyyzz[k] = -g_z_0_zzzz_yyyyzz[k] * cd_x[k] + g_z_0_zzzz_xyyyyzz[k];

                g_z_0_xzzzz_yyyzzz[k] = -g_z_0_zzzz_yyyzzz[k] * cd_x[k] + g_z_0_zzzz_xyyyzzz[k];

                g_z_0_xzzzz_yyzzzz[k] = -g_z_0_zzzz_yyzzzz[k] * cd_x[k] + g_z_0_zzzz_xyyzzzz[k];

                g_z_0_xzzzz_yzzzzz[k] = -g_z_0_zzzz_yzzzzz[k] * cd_x[k] + g_z_0_zzzz_xyzzzzz[k];

                g_z_0_xzzzz_zzzzzz[k] = -g_z_0_zzzz_zzzzzz[k] * cd_x[k] + g_z_0_zzzz_xzzzzzz[k];
            }

            /// Set up 420-448 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 420);

            auto g_z_0_yyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 421);

            auto g_z_0_yyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 422);

            auto g_z_0_yyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 423);

            auto g_z_0_yyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 424);

            auto g_z_0_yyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 425);

            auto g_z_0_yyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 426);

            auto g_z_0_yyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 427);

            auto g_z_0_yyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 428);

            auto g_z_0_yyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 429);

            auto g_z_0_yyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 430);

            auto g_z_0_yyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 431);

            auto g_z_0_yyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 432);

            auto g_z_0_yyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 433);

            auto g_z_0_yyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 434);

            auto g_z_0_yyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 435);

            auto g_z_0_yyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 436);

            auto g_z_0_yyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 437);

            auto g_z_0_yyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 438);

            auto g_z_0_yyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 439);

            auto g_z_0_yyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 440);

            auto g_z_0_yyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 441);

            auto g_z_0_yyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 442);

            auto g_z_0_yyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 443);

            auto g_z_0_yyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 444);

            auto g_z_0_yyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 445);

            auto g_z_0_yyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 446);

            auto g_z_0_yyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 447);

            #pragma omp simd aligned(cd_y, g_z_0_yyyy_xxxxxx, g_z_0_yyyy_xxxxxxy, g_z_0_yyyy_xxxxxy, g_z_0_yyyy_xxxxxyy, g_z_0_yyyy_xxxxxyz, g_z_0_yyyy_xxxxxz, g_z_0_yyyy_xxxxyy, g_z_0_yyyy_xxxxyyy, g_z_0_yyyy_xxxxyyz, g_z_0_yyyy_xxxxyz, g_z_0_yyyy_xxxxyzz, g_z_0_yyyy_xxxxzz, g_z_0_yyyy_xxxyyy, g_z_0_yyyy_xxxyyyy, g_z_0_yyyy_xxxyyyz, g_z_0_yyyy_xxxyyz, g_z_0_yyyy_xxxyyzz, g_z_0_yyyy_xxxyzz, g_z_0_yyyy_xxxyzzz, g_z_0_yyyy_xxxzzz, g_z_0_yyyy_xxyyyy, g_z_0_yyyy_xxyyyyy, g_z_0_yyyy_xxyyyyz, g_z_0_yyyy_xxyyyz, g_z_0_yyyy_xxyyyzz, g_z_0_yyyy_xxyyzz, g_z_0_yyyy_xxyyzzz, g_z_0_yyyy_xxyzzz, g_z_0_yyyy_xxyzzzz, g_z_0_yyyy_xxzzzz, g_z_0_yyyy_xyyyyy, g_z_0_yyyy_xyyyyyy, g_z_0_yyyy_xyyyyyz, g_z_0_yyyy_xyyyyz, g_z_0_yyyy_xyyyyzz, g_z_0_yyyy_xyyyzz, g_z_0_yyyy_xyyyzzz, g_z_0_yyyy_xyyzzz, g_z_0_yyyy_xyyzzzz, g_z_0_yyyy_xyzzzz, g_z_0_yyyy_xyzzzzz, g_z_0_yyyy_xzzzzz, g_z_0_yyyy_yyyyyy, g_z_0_yyyy_yyyyyyy, g_z_0_yyyy_yyyyyyz, g_z_0_yyyy_yyyyyz, g_z_0_yyyy_yyyyyzz, g_z_0_yyyy_yyyyzz, g_z_0_yyyy_yyyyzzz, g_z_0_yyyy_yyyzzz, g_z_0_yyyy_yyyzzzz, g_z_0_yyyy_yyzzzz, g_z_0_yyyy_yyzzzzz, g_z_0_yyyy_yzzzzz, g_z_0_yyyy_yzzzzzz, g_z_0_yyyy_zzzzzz, g_z_0_yyyyy_xxxxxx, g_z_0_yyyyy_xxxxxy, g_z_0_yyyyy_xxxxxz, g_z_0_yyyyy_xxxxyy, g_z_0_yyyyy_xxxxyz, g_z_0_yyyyy_xxxxzz, g_z_0_yyyyy_xxxyyy, g_z_0_yyyyy_xxxyyz, g_z_0_yyyyy_xxxyzz, g_z_0_yyyyy_xxxzzz, g_z_0_yyyyy_xxyyyy, g_z_0_yyyyy_xxyyyz, g_z_0_yyyyy_xxyyzz, g_z_0_yyyyy_xxyzzz, g_z_0_yyyyy_xxzzzz, g_z_0_yyyyy_xyyyyy, g_z_0_yyyyy_xyyyyz, g_z_0_yyyyy_xyyyzz, g_z_0_yyyyy_xyyzzz, g_z_0_yyyyy_xyzzzz, g_z_0_yyyyy_xzzzzz, g_z_0_yyyyy_yyyyyy, g_z_0_yyyyy_yyyyyz, g_z_0_yyyyy_yyyyzz, g_z_0_yyyyy_yyyzzz, g_z_0_yyyyy_yyzzzz, g_z_0_yyyyy_yzzzzz, g_z_0_yyyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyy_xxxxxx[k] = -g_z_0_yyyy_xxxxxx[k] * cd_y[k] + g_z_0_yyyy_xxxxxxy[k];

                g_z_0_yyyyy_xxxxxy[k] = -g_z_0_yyyy_xxxxxy[k] * cd_y[k] + g_z_0_yyyy_xxxxxyy[k];

                g_z_0_yyyyy_xxxxxz[k] = -g_z_0_yyyy_xxxxxz[k] * cd_y[k] + g_z_0_yyyy_xxxxxyz[k];

                g_z_0_yyyyy_xxxxyy[k] = -g_z_0_yyyy_xxxxyy[k] * cd_y[k] + g_z_0_yyyy_xxxxyyy[k];

                g_z_0_yyyyy_xxxxyz[k] = -g_z_0_yyyy_xxxxyz[k] * cd_y[k] + g_z_0_yyyy_xxxxyyz[k];

                g_z_0_yyyyy_xxxxzz[k] = -g_z_0_yyyy_xxxxzz[k] * cd_y[k] + g_z_0_yyyy_xxxxyzz[k];

                g_z_0_yyyyy_xxxyyy[k] = -g_z_0_yyyy_xxxyyy[k] * cd_y[k] + g_z_0_yyyy_xxxyyyy[k];

                g_z_0_yyyyy_xxxyyz[k] = -g_z_0_yyyy_xxxyyz[k] * cd_y[k] + g_z_0_yyyy_xxxyyyz[k];

                g_z_0_yyyyy_xxxyzz[k] = -g_z_0_yyyy_xxxyzz[k] * cd_y[k] + g_z_0_yyyy_xxxyyzz[k];

                g_z_0_yyyyy_xxxzzz[k] = -g_z_0_yyyy_xxxzzz[k] * cd_y[k] + g_z_0_yyyy_xxxyzzz[k];

                g_z_0_yyyyy_xxyyyy[k] = -g_z_0_yyyy_xxyyyy[k] * cd_y[k] + g_z_0_yyyy_xxyyyyy[k];

                g_z_0_yyyyy_xxyyyz[k] = -g_z_0_yyyy_xxyyyz[k] * cd_y[k] + g_z_0_yyyy_xxyyyyz[k];

                g_z_0_yyyyy_xxyyzz[k] = -g_z_0_yyyy_xxyyzz[k] * cd_y[k] + g_z_0_yyyy_xxyyyzz[k];

                g_z_0_yyyyy_xxyzzz[k] = -g_z_0_yyyy_xxyzzz[k] * cd_y[k] + g_z_0_yyyy_xxyyzzz[k];

                g_z_0_yyyyy_xxzzzz[k] = -g_z_0_yyyy_xxzzzz[k] * cd_y[k] + g_z_0_yyyy_xxyzzzz[k];

                g_z_0_yyyyy_xyyyyy[k] = -g_z_0_yyyy_xyyyyy[k] * cd_y[k] + g_z_0_yyyy_xyyyyyy[k];

                g_z_0_yyyyy_xyyyyz[k] = -g_z_0_yyyy_xyyyyz[k] * cd_y[k] + g_z_0_yyyy_xyyyyyz[k];

                g_z_0_yyyyy_xyyyzz[k] = -g_z_0_yyyy_xyyyzz[k] * cd_y[k] + g_z_0_yyyy_xyyyyzz[k];

                g_z_0_yyyyy_xyyzzz[k] = -g_z_0_yyyy_xyyzzz[k] * cd_y[k] + g_z_0_yyyy_xyyyzzz[k];

                g_z_0_yyyyy_xyzzzz[k] = -g_z_0_yyyy_xyzzzz[k] * cd_y[k] + g_z_0_yyyy_xyyzzzz[k];

                g_z_0_yyyyy_xzzzzz[k] = -g_z_0_yyyy_xzzzzz[k] * cd_y[k] + g_z_0_yyyy_xyzzzzz[k];

                g_z_0_yyyyy_yyyyyy[k] = -g_z_0_yyyy_yyyyyy[k] * cd_y[k] + g_z_0_yyyy_yyyyyyy[k];

                g_z_0_yyyyy_yyyyyz[k] = -g_z_0_yyyy_yyyyyz[k] * cd_y[k] + g_z_0_yyyy_yyyyyyz[k];

                g_z_0_yyyyy_yyyyzz[k] = -g_z_0_yyyy_yyyyzz[k] * cd_y[k] + g_z_0_yyyy_yyyyyzz[k];

                g_z_0_yyyyy_yyyzzz[k] = -g_z_0_yyyy_yyyzzz[k] * cd_y[k] + g_z_0_yyyy_yyyyzzz[k];

                g_z_0_yyyyy_yyzzzz[k] = -g_z_0_yyyy_yyzzzz[k] * cd_y[k] + g_z_0_yyyy_yyyzzzz[k];

                g_z_0_yyyyy_yzzzzz[k] = -g_z_0_yyyy_yzzzzz[k] * cd_y[k] + g_z_0_yyyy_yyzzzzz[k];

                g_z_0_yyyyy_zzzzzz[k] = -g_z_0_yyyy_zzzzzz[k] * cd_y[k] + g_z_0_yyyy_yzzzzzz[k];
            }

            /// Set up 448-476 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 448);

            auto g_z_0_yyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 449);

            auto g_z_0_yyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 450);

            auto g_z_0_yyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 451);

            auto g_z_0_yyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 452);

            auto g_z_0_yyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 453);

            auto g_z_0_yyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 454);

            auto g_z_0_yyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 455);

            auto g_z_0_yyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 456);

            auto g_z_0_yyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 457);

            auto g_z_0_yyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 458);

            auto g_z_0_yyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 459);

            auto g_z_0_yyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 460);

            auto g_z_0_yyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 461);

            auto g_z_0_yyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 462);

            auto g_z_0_yyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 463);

            auto g_z_0_yyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 464);

            auto g_z_0_yyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 465);

            auto g_z_0_yyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 466);

            auto g_z_0_yyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 467);

            auto g_z_0_yyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 468);

            auto g_z_0_yyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 469);

            auto g_z_0_yyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 470);

            auto g_z_0_yyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 471);

            auto g_z_0_yyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 472);

            auto g_z_0_yyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 473);

            auto g_z_0_yyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 474);

            auto g_z_0_yyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 475);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyz_xxxxxx, g_z_0_yyyyz_xxxxxy, g_z_0_yyyyz_xxxxxz, g_z_0_yyyyz_xxxxyy, g_z_0_yyyyz_xxxxyz, g_z_0_yyyyz_xxxxzz, g_z_0_yyyyz_xxxyyy, g_z_0_yyyyz_xxxyyz, g_z_0_yyyyz_xxxyzz, g_z_0_yyyyz_xxxzzz, g_z_0_yyyyz_xxyyyy, g_z_0_yyyyz_xxyyyz, g_z_0_yyyyz_xxyyzz, g_z_0_yyyyz_xxyzzz, g_z_0_yyyyz_xxzzzz, g_z_0_yyyyz_xyyyyy, g_z_0_yyyyz_xyyyyz, g_z_0_yyyyz_xyyyzz, g_z_0_yyyyz_xyyzzz, g_z_0_yyyyz_xyzzzz, g_z_0_yyyyz_xzzzzz, g_z_0_yyyyz_yyyyyy, g_z_0_yyyyz_yyyyyz, g_z_0_yyyyz_yyyyzz, g_z_0_yyyyz_yyyzzz, g_z_0_yyyyz_yyzzzz, g_z_0_yyyyz_yzzzzz, g_z_0_yyyyz_zzzzzz, g_z_0_yyyz_xxxxxx, g_z_0_yyyz_xxxxxxy, g_z_0_yyyz_xxxxxy, g_z_0_yyyz_xxxxxyy, g_z_0_yyyz_xxxxxyz, g_z_0_yyyz_xxxxxz, g_z_0_yyyz_xxxxyy, g_z_0_yyyz_xxxxyyy, g_z_0_yyyz_xxxxyyz, g_z_0_yyyz_xxxxyz, g_z_0_yyyz_xxxxyzz, g_z_0_yyyz_xxxxzz, g_z_0_yyyz_xxxyyy, g_z_0_yyyz_xxxyyyy, g_z_0_yyyz_xxxyyyz, g_z_0_yyyz_xxxyyz, g_z_0_yyyz_xxxyyzz, g_z_0_yyyz_xxxyzz, g_z_0_yyyz_xxxyzzz, g_z_0_yyyz_xxxzzz, g_z_0_yyyz_xxyyyy, g_z_0_yyyz_xxyyyyy, g_z_0_yyyz_xxyyyyz, g_z_0_yyyz_xxyyyz, g_z_0_yyyz_xxyyyzz, g_z_0_yyyz_xxyyzz, g_z_0_yyyz_xxyyzzz, g_z_0_yyyz_xxyzzz, g_z_0_yyyz_xxyzzzz, g_z_0_yyyz_xxzzzz, g_z_0_yyyz_xyyyyy, g_z_0_yyyz_xyyyyyy, g_z_0_yyyz_xyyyyyz, g_z_0_yyyz_xyyyyz, g_z_0_yyyz_xyyyyzz, g_z_0_yyyz_xyyyzz, g_z_0_yyyz_xyyyzzz, g_z_0_yyyz_xyyzzz, g_z_0_yyyz_xyyzzzz, g_z_0_yyyz_xyzzzz, g_z_0_yyyz_xyzzzzz, g_z_0_yyyz_xzzzzz, g_z_0_yyyz_yyyyyy, g_z_0_yyyz_yyyyyyy, g_z_0_yyyz_yyyyyyz, g_z_0_yyyz_yyyyyz, g_z_0_yyyz_yyyyyzz, g_z_0_yyyz_yyyyzz, g_z_0_yyyz_yyyyzzz, g_z_0_yyyz_yyyzzz, g_z_0_yyyz_yyyzzzz, g_z_0_yyyz_yyzzzz, g_z_0_yyyz_yyzzzzz, g_z_0_yyyz_yzzzzz, g_z_0_yyyz_yzzzzzz, g_z_0_yyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyz_xxxxxx[k] = -g_z_0_yyyz_xxxxxx[k] * cd_y[k] + g_z_0_yyyz_xxxxxxy[k];

                g_z_0_yyyyz_xxxxxy[k] = -g_z_0_yyyz_xxxxxy[k] * cd_y[k] + g_z_0_yyyz_xxxxxyy[k];

                g_z_0_yyyyz_xxxxxz[k] = -g_z_0_yyyz_xxxxxz[k] * cd_y[k] + g_z_0_yyyz_xxxxxyz[k];

                g_z_0_yyyyz_xxxxyy[k] = -g_z_0_yyyz_xxxxyy[k] * cd_y[k] + g_z_0_yyyz_xxxxyyy[k];

                g_z_0_yyyyz_xxxxyz[k] = -g_z_0_yyyz_xxxxyz[k] * cd_y[k] + g_z_0_yyyz_xxxxyyz[k];

                g_z_0_yyyyz_xxxxzz[k] = -g_z_0_yyyz_xxxxzz[k] * cd_y[k] + g_z_0_yyyz_xxxxyzz[k];

                g_z_0_yyyyz_xxxyyy[k] = -g_z_0_yyyz_xxxyyy[k] * cd_y[k] + g_z_0_yyyz_xxxyyyy[k];

                g_z_0_yyyyz_xxxyyz[k] = -g_z_0_yyyz_xxxyyz[k] * cd_y[k] + g_z_0_yyyz_xxxyyyz[k];

                g_z_0_yyyyz_xxxyzz[k] = -g_z_0_yyyz_xxxyzz[k] * cd_y[k] + g_z_0_yyyz_xxxyyzz[k];

                g_z_0_yyyyz_xxxzzz[k] = -g_z_0_yyyz_xxxzzz[k] * cd_y[k] + g_z_0_yyyz_xxxyzzz[k];

                g_z_0_yyyyz_xxyyyy[k] = -g_z_0_yyyz_xxyyyy[k] * cd_y[k] + g_z_0_yyyz_xxyyyyy[k];

                g_z_0_yyyyz_xxyyyz[k] = -g_z_0_yyyz_xxyyyz[k] * cd_y[k] + g_z_0_yyyz_xxyyyyz[k];

                g_z_0_yyyyz_xxyyzz[k] = -g_z_0_yyyz_xxyyzz[k] * cd_y[k] + g_z_0_yyyz_xxyyyzz[k];

                g_z_0_yyyyz_xxyzzz[k] = -g_z_0_yyyz_xxyzzz[k] * cd_y[k] + g_z_0_yyyz_xxyyzzz[k];

                g_z_0_yyyyz_xxzzzz[k] = -g_z_0_yyyz_xxzzzz[k] * cd_y[k] + g_z_0_yyyz_xxyzzzz[k];

                g_z_0_yyyyz_xyyyyy[k] = -g_z_0_yyyz_xyyyyy[k] * cd_y[k] + g_z_0_yyyz_xyyyyyy[k];

                g_z_0_yyyyz_xyyyyz[k] = -g_z_0_yyyz_xyyyyz[k] * cd_y[k] + g_z_0_yyyz_xyyyyyz[k];

                g_z_0_yyyyz_xyyyzz[k] = -g_z_0_yyyz_xyyyzz[k] * cd_y[k] + g_z_0_yyyz_xyyyyzz[k];

                g_z_0_yyyyz_xyyzzz[k] = -g_z_0_yyyz_xyyzzz[k] * cd_y[k] + g_z_0_yyyz_xyyyzzz[k];

                g_z_0_yyyyz_xyzzzz[k] = -g_z_0_yyyz_xyzzzz[k] * cd_y[k] + g_z_0_yyyz_xyyzzzz[k];

                g_z_0_yyyyz_xzzzzz[k] = -g_z_0_yyyz_xzzzzz[k] * cd_y[k] + g_z_0_yyyz_xyzzzzz[k];

                g_z_0_yyyyz_yyyyyy[k] = -g_z_0_yyyz_yyyyyy[k] * cd_y[k] + g_z_0_yyyz_yyyyyyy[k];

                g_z_0_yyyyz_yyyyyz[k] = -g_z_0_yyyz_yyyyyz[k] * cd_y[k] + g_z_0_yyyz_yyyyyyz[k];

                g_z_0_yyyyz_yyyyzz[k] = -g_z_0_yyyz_yyyyzz[k] * cd_y[k] + g_z_0_yyyz_yyyyyzz[k];

                g_z_0_yyyyz_yyyzzz[k] = -g_z_0_yyyz_yyyzzz[k] * cd_y[k] + g_z_0_yyyz_yyyyzzz[k];

                g_z_0_yyyyz_yyzzzz[k] = -g_z_0_yyyz_yyzzzz[k] * cd_y[k] + g_z_0_yyyz_yyyzzzz[k];

                g_z_0_yyyyz_yzzzzz[k] = -g_z_0_yyyz_yzzzzz[k] * cd_y[k] + g_z_0_yyyz_yyzzzzz[k];

                g_z_0_yyyyz_zzzzzz[k] = -g_z_0_yyyz_zzzzzz[k] * cd_y[k] + g_z_0_yyyz_yzzzzzz[k];
            }

            /// Set up 476-504 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 476);

            auto g_z_0_yyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 477);

            auto g_z_0_yyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 478);

            auto g_z_0_yyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 479);

            auto g_z_0_yyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 480);

            auto g_z_0_yyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 481);

            auto g_z_0_yyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 482);

            auto g_z_0_yyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 483);

            auto g_z_0_yyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 484);

            auto g_z_0_yyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 485);

            auto g_z_0_yyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 486);

            auto g_z_0_yyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 487);

            auto g_z_0_yyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 488);

            auto g_z_0_yyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 489);

            auto g_z_0_yyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 490);

            auto g_z_0_yyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 491);

            auto g_z_0_yyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 492);

            auto g_z_0_yyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 493);

            auto g_z_0_yyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 494);

            auto g_z_0_yyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 495);

            auto g_z_0_yyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 496);

            auto g_z_0_yyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 497);

            auto g_z_0_yyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 498);

            auto g_z_0_yyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 499);

            auto g_z_0_yyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 500);

            auto g_z_0_yyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 501);

            auto g_z_0_yyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 502);

            auto g_z_0_yyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 503);

            #pragma omp simd aligned(cd_y, g_z_0_yyyzz_xxxxxx, g_z_0_yyyzz_xxxxxy, g_z_0_yyyzz_xxxxxz, g_z_0_yyyzz_xxxxyy, g_z_0_yyyzz_xxxxyz, g_z_0_yyyzz_xxxxzz, g_z_0_yyyzz_xxxyyy, g_z_0_yyyzz_xxxyyz, g_z_0_yyyzz_xxxyzz, g_z_0_yyyzz_xxxzzz, g_z_0_yyyzz_xxyyyy, g_z_0_yyyzz_xxyyyz, g_z_0_yyyzz_xxyyzz, g_z_0_yyyzz_xxyzzz, g_z_0_yyyzz_xxzzzz, g_z_0_yyyzz_xyyyyy, g_z_0_yyyzz_xyyyyz, g_z_0_yyyzz_xyyyzz, g_z_0_yyyzz_xyyzzz, g_z_0_yyyzz_xyzzzz, g_z_0_yyyzz_xzzzzz, g_z_0_yyyzz_yyyyyy, g_z_0_yyyzz_yyyyyz, g_z_0_yyyzz_yyyyzz, g_z_0_yyyzz_yyyzzz, g_z_0_yyyzz_yyzzzz, g_z_0_yyyzz_yzzzzz, g_z_0_yyyzz_zzzzzz, g_z_0_yyzz_xxxxxx, g_z_0_yyzz_xxxxxxy, g_z_0_yyzz_xxxxxy, g_z_0_yyzz_xxxxxyy, g_z_0_yyzz_xxxxxyz, g_z_0_yyzz_xxxxxz, g_z_0_yyzz_xxxxyy, g_z_0_yyzz_xxxxyyy, g_z_0_yyzz_xxxxyyz, g_z_0_yyzz_xxxxyz, g_z_0_yyzz_xxxxyzz, g_z_0_yyzz_xxxxzz, g_z_0_yyzz_xxxyyy, g_z_0_yyzz_xxxyyyy, g_z_0_yyzz_xxxyyyz, g_z_0_yyzz_xxxyyz, g_z_0_yyzz_xxxyyzz, g_z_0_yyzz_xxxyzz, g_z_0_yyzz_xxxyzzz, g_z_0_yyzz_xxxzzz, g_z_0_yyzz_xxyyyy, g_z_0_yyzz_xxyyyyy, g_z_0_yyzz_xxyyyyz, g_z_0_yyzz_xxyyyz, g_z_0_yyzz_xxyyyzz, g_z_0_yyzz_xxyyzz, g_z_0_yyzz_xxyyzzz, g_z_0_yyzz_xxyzzz, g_z_0_yyzz_xxyzzzz, g_z_0_yyzz_xxzzzz, g_z_0_yyzz_xyyyyy, g_z_0_yyzz_xyyyyyy, g_z_0_yyzz_xyyyyyz, g_z_0_yyzz_xyyyyz, g_z_0_yyzz_xyyyyzz, g_z_0_yyzz_xyyyzz, g_z_0_yyzz_xyyyzzz, g_z_0_yyzz_xyyzzz, g_z_0_yyzz_xyyzzzz, g_z_0_yyzz_xyzzzz, g_z_0_yyzz_xyzzzzz, g_z_0_yyzz_xzzzzz, g_z_0_yyzz_yyyyyy, g_z_0_yyzz_yyyyyyy, g_z_0_yyzz_yyyyyyz, g_z_0_yyzz_yyyyyz, g_z_0_yyzz_yyyyyzz, g_z_0_yyzz_yyyyzz, g_z_0_yyzz_yyyyzzz, g_z_0_yyzz_yyyzzz, g_z_0_yyzz_yyyzzzz, g_z_0_yyzz_yyzzzz, g_z_0_yyzz_yyzzzzz, g_z_0_yyzz_yzzzzz, g_z_0_yyzz_yzzzzzz, g_z_0_yyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzz_xxxxxx[k] = -g_z_0_yyzz_xxxxxx[k] * cd_y[k] + g_z_0_yyzz_xxxxxxy[k];

                g_z_0_yyyzz_xxxxxy[k] = -g_z_0_yyzz_xxxxxy[k] * cd_y[k] + g_z_0_yyzz_xxxxxyy[k];

                g_z_0_yyyzz_xxxxxz[k] = -g_z_0_yyzz_xxxxxz[k] * cd_y[k] + g_z_0_yyzz_xxxxxyz[k];

                g_z_0_yyyzz_xxxxyy[k] = -g_z_0_yyzz_xxxxyy[k] * cd_y[k] + g_z_0_yyzz_xxxxyyy[k];

                g_z_0_yyyzz_xxxxyz[k] = -g_z_0_yyzz_xxxxyz[k] * cd_y[k] + g_z_0_yyzz_xxxxyyz[k];

                g_z_0_yyyzz_xxxxzz[k] = -g_z_0_yyzz_xxxxzz[k] * cd_y[k] + g_z_0_yyzz_xxxxyzz[k];

                g_z_0_yyyzz_xxxyyy[k] = -g_z_0_yyzz_xxxyyy[k] * cd_y[k] + g_z_0_yyzz_xxxyyyy[k];

                g_z_0_yyyzz_xxxyyz[k] = -g_z_0_yyzz_xxxyyz[k] * cd_y[k] + g_z_0_yyzz_xxxyyyz[k];

                g_z_0_yyyzz_xxxyzz[k] = -g_z_0_yyzz_xxxyzz[k] * cd_y[k] + g_z_0_yyzz_xxxyyzz[k];

                g_z_0_yyyzz_xxxzzz[k] = -g_z_0_yyzz_xxxzzz[k] * cd_y[k] + g_z_0_yyzz_xxxyzzz[k];

                g_z_0_yyyzz_xxyyyy[k] = -g_z_0_yyzz_xxyyyy[k] * cd_y[k] + g_z_0_yyzz_xxyyyyy[k];

                g_z_0_yyyzz_xxyyyz[k] = -g_z_0_yyzz_xxyyyz[k] * cd_y[k] + g_z_0_yyzz_xxyyyyz[k];

                g_z_0_yyyzz_xxyyzz[k] = -g_z_0_yyzz_xxyyzz[k] * cd_y[k] + g_z_0_yyzz_xxyyyzz[k];

                g_z_0_yyyzz_xxyzzz[k] = -g_z_0_yyzz_xxyzzz[k] * cd_y[k] + g_z_0_yyzz_xxyyzzz[k];

                g_z_0_yyyzz_xxzzzz[k] = -g_z_0_yyzz_xxzzzz[k] * cd_y[k] + g_z_0_yyzz_xxyzzzz[k];

                g_z_0_yyyzz_xyyyyy[k] = -g_z_0_yyzz_xyyyyy[k] * cd_y[k] + g_z_0_yyzz_xyyyyyy[k];

                g_z_0_yyyzz_xyyyyz[k] = -g_z_0_yyzz_xyyyyz[k] * cd_y[k] + g_z_0_yyzz_xyyyyyz[k];

                g_z_0_yyyzz_xyyyzz[k] = -g_z_0_yyzz_xyyyzz[k] * cd_y[k] + g_z_0_yyzz_xyyyyzz[k];

                g_z_0_yyyzz_xyyzzz[k] = -g_z_0_yyzz_xyyzzz[k] * cd_y[k] + g_z_0_yyzz_xyyyzzz[k];

                g_z_0_yyyzz_xyzzzz[k] = -g_z_0_yyzz_xyzzzz[k] * cd_y[k] + g_z_0_yyzz_xyyzzzz[k];

                g_z_0_yyyzz_xzzzzz[k] = -g_z_0_yyzz_xzzzzz[k] * cd_y[k] + g_z_0_yyzz_xyzzzzz[k];

                g_z_0_yyyzz_yyyyyy[k] = -g_z_0_yyzz_yyyyyy[k] * cd_y[k] + g_z_0_yyzz_yyyyyyy[k];

                g_z_0_yyyzz_yyyyyz[k] = -g_z_0_yyzz_yyyyyz[k] * cd_y[k] + g_z_0_yyzz_yyyyyyz[k];

                g_z_0_yyyzz_yyyyzz[k] = -g_z_0_yyzz_yyyyzz[k] * cd_y[k] + g_z_0_yyzz_yyyyyzz[k];

                g_z_0_yyyzz_yyyzzz[k] = -g_z_0_yyzz_yyyzzz[k] * cd_y[k] + g_z_0_yyzz_yyyyzzz[k];

                g_z_0_yyyzz_yyzzzz[k] = -g_z_0_yyzz_yyzzzz[k] * cd_y[k] + g_z_0_yyzz_yyyzzzz[k];

                g_z_0_yyyzz_yzzzzz[k] = -g_z_0_yyzz_yzzzzz[k] * cd_y[k] + g_z_0_yyzz_yyzzzzz[k];

                g_z_0_yyyzz_zzzzzz[k] = -g_z_0_yyzz_zzzzzz[k] * cd_y[k] + g_z_0_yyzz_yzzzzzz[k];
            }

            /// Set up 504-532 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 504);

            auto g_z_0_yyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 505);

            auto g_z_0_yyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 506);

            auto g_z_0_yyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 507);

            auto g_z_0_yyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 508);

            auto g_z_0_yyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 509);

            auto g_z_0_yyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 510);

            auto g_z_0_yyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 511);

            auto g_z_0_yyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 512);

            auto g_z_0_yyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 513);

            auto g_z_0_yyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 514);

            auto g_z_0_yyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 515);

            auto g_z_0_yyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 516);

            auto g_z_0_yyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 517);

            auto g_z_0_yyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 518);

            auto g_z_0_yyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 519);

            auto g_z_0_yyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 520);

            auto g_z_0_yyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 521);

            auto g_z_0_yyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 522);

            auto g_z_0_yyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 523);

            auto g_z_0_yyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 524);

            auto g_z_0_yyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 525);

            auto g_z_0_yyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 526);

            auto g_z_0_yyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 527);

            auto g_z_0_yyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 528);

            auto g_z_0_yyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 529);

            auto g_z_0_yyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 530);

            auto g_z_0_yyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 531);

            #pragma omp simd aligned(cd_y, g_z_0_yyzzz_xxxxxx, g_z_0_yyzzz_xxxxxy, g_z_0_yyzzz_xxxxxz, g_z_0_yyzzz_xxxxyy, g_z_0_yyzzz_xxxxyz, g_z_0_yyzzz_xxxxzz, g_z_0_yyzzz_xxxyyy, g_z_0_yyzzz_xxxyyz, g_z_0_yyzzz_xxxyzz, g_z_0_yyzzz_xxxzzz, g_z_0_yyzzz_xxyyyy, g_z_0_yyzzz_xxyyyz, g_z_0_yyzzz_xxyyzz, g_z_0_yyzzz_xxyzzz, g_z_0_yyzzz_xxzzzz, g_z_0_yyzzz_xyyyyy, g_z_0_yyzzz_xyyyyz, g_z_0_yyzzz_xyyyzz, g_z_0_yyzzz_xyyzzz, g_z_0_yyzzz_xyzzzz, g_z_0_yyzzz_xzzzzz, g_z_0_yyzzz_yyyyyy, g_z_0_yyzzz_yyyyyz, g_z_0_yyzzz_yyyyzz, g_z_0_yyzzz_yyyzzz, g_z_0_yyzzz_yyzzzz, g_z_0_yyzzz_yzzzzz, g_z_0_yyzzz_zzzzzz, g_z_0_yzzz_xxxxxx, g_z_0_yzzz_xxxxxxy, g_z_0_yzzz_xxxxxy, g_z_0_yzzz_xxxxxyy, g_z_0_yzzz_xxxxxyz, g_z_0_yzzz_xxxxxz, g_z_0_yzzz_xxxxyy, g_z_0_yzzz_xxxxyyy, g_z_0_yzzz_xxxxyyz, g_z_0_yzzz_xxxxyz, g_z_0_yzzz_xxxxyzz, g_z_0_yzzz_xxxxzz, g_z_0_yzzz_xxxyyy, g_z_0_yzzz_xxxyyyy, g_z_0_yzzz_xxxyyyz, g_z_0_yzzz_xxxyyz, g_z_0_yzzz_xxxyyzz, g_z_0_yzzz_xxxyzz, g_z_0_yzzz_xxxyzzz, g_z_0_yzzz_xxxzzz, g_z_0_yzzz_xxyyyy, g_z_0_yzzz_xxyyyyy, g_z_0_yzzz_xxyyyyz, g_z_0_yzzz_xxyyyz, g_z_0_yzzz_xxyyyzz, g_z_0_yzzz_xxyyzz, g_z_0_yzzz_xxyyzzz, g_z_0_yzzz_xxyzzz, g_z_0_yzzz_xxyzzzz, g_z_0_yzzz_xxzzzz, g_z_0_yzzz_xyyyyy, g_z_0_yzzz_xyyyyyy, g_z_0_yzzz_xyyyyyz, g_z_0_yzzz_xyyyyz, g_z_0_yzzz_xyyyyzz, g_z_0_yzzz_xyyyzz, g_z_0_yzzz_xyyyzzz, g_z_0_yzzz_xyyzzz, g_z_0_yzzz_xyyzzzz, g_z_0_yzzz_xyzzzz, g_z_0_yzzz_xyzzzzz, g_z_0_yzzz_xzzzzz, g_z_0_yzzz_yyyyyy, g_z_0_yzzz_yyyyyyy, g_z_0_yzzz_yyyyyyz, g_z_0_yzzz_yyyyyz, g_z_0_yzzz_yyyyyzz, g_z_0_yzzz_yyyyzz, g_z_0_yzzz_yyyyzzz, g_z_0_yzzz_yyyzzz, g_z_0_yzzz_yyyzzzz, g_z_0_yzzz_yyzzzz, g_z_0_yzzz_yyzzzzz, g_z_0_yzzz_yzzzzz, g_z_0_yzzz_yzzzzzz, g_z_0_yzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzz_xxxxxx[k] = -g_z_0_yzzz_xxxxxx[k] * cd_y[k] + g_z_0_yzzz_xxxxxxy[k];

                g_z_0_yyzzz_xxxxxy[k] = -g_z_0_yzzz_xxxxxy[k] * cd_y[k] + g_z_0_yzzz_xxxxxyy[k];

                g_z_0_yyzzz_xxxxxz[k] = -g_z_0_yzzz_xxxxxz[k] * cd_y[k] + g_z_0_yzzz_xxxxxyz[k];

                g_z_0_yyzzz_xxxxyy[k] = -g_z_0_yzzz_xxxxyy[k] * cd_y[k] + g_z_0_yzzz_xxxxyyy[k];

                g_z_0_yyzzz_xxxxyz[k] = -g_z_0_yzzz_xxxxyz[k] * cd_y[k] + g_z_0_yzzz_xxxxyyz[k];

                g_z_0_yyzzz_xxxxzz[k] = -g_z_0_yzzz_xxxxzz[k] * cd_y[k] + g_z_0_yzzz_xxxxyzz[k];

                g_z_0_yyzzz_xxxyyy[k] = -g_z_0_yzzz_xxxyyy[k] * cd_y[k] + g_z_0_yzzz_xxxyyyy[k];

                g_z_0_yyzzz_xxxyyz[k] = -g_z_0_yzzz_xxxyyz[k] * cd_y[k] + g_z_0_yzzz_xxxyyyz[k];

                g_z_0_yyzzz_xxxyzz[k] = -g_z_0_yzzz_xxxyzz[k] * cd_y[k] + g_z_0_yzzz_xxxyyzz[k];

                g_z_0_yyzzz_xxxzzz[k] = -g_z_0_yzzz_xxxzzz[k] * cd_y[k] + g_z_0_yzzz_xxxyzzz[k];

                g_z_0_yyzzz_xxyyyy[k] = -g_z_0_yzzz_xxyyyy[k] * cd_y[k] + g_z_0_yzzz_xxyyyyy[k];

                g_z_0_yyzzz_xxyyyz[k] = -g_z_0_yzzz_xxyyyz[k] * cd_y[k] + g_z_0_yzzz_xxyyyyz[k];

                g_z_0_yyzzz_xxyyzz[k] = -g_z_0_yzzz_xxyyzz[k] * cd_y[k] + g_z_0_yzzz_xxyyyzz[k];

                g_z_0_yyzzz_xxyzzz[k] = -g_z_0_yzzz_xxyzzz[k] * cd_y[k] + g_z_0_yzzz_xxyyzzz[k];

                g_z_0_yyzzz_xxzzzz[k] = -g_z_0_yzzz_xxzzzz[k] * cd_y[k] + g_z_0_yzzz_xxyzzzz[k];

                g_z_0_yyzzz_xyyyyy[k] = -g_z_0_yzzz_xyyyyy[k] * cd_y[k] + g_z_0_yzzz_xyyyyyy[k];

                g_z_0_yyzzz_xyyyyz[k] = -g_z_0_yzzz_xyyyyz[k] * cd_y[k] + g_z_0_yzzz_xyyyyyz[k];

                g_z_0_yyzzz_xyyyzz[k] = -g_z_0_yzzz_xyyyzz[k] * cd_y[k] + g_z_0_yzzz_xyyyyzz[k];

                g_z_0_yyzzz_xyyzzz[k] = -g_z_0_yzzz_xyyzzz[k] * cd_y[k] + g_z_0_yzzz_xyyyzzz[k];

                g_z_0_yyzzz_xyzzzz[k] = -g_z_0_yzzz_xyzzzz[k] * cd_y[k] + g_z_0_yzzz_xyyzzzz[k];

                g_z_0_yyzzz_xzzzzz[k] = -g_z_0_yzzz_xzzzzz[k] * cd_y[k] + g_z_0_yzzz_xyzzzzz[k];

                g_z_0_yyzzz_yyyyyy[k] = -g_z_0_yzzz_yyyyyy[k] * cd_y[k] + g_z_0_yzzz_yyyyyyy[k];

                g_z_0_yyzzz_yyyyyz[k] = -g_z_0_yzzz_yyyyyz[k] * cd_y[k] + g_z_0_yzzz_yyyyyyz[k];

                g_z_0_yyzzz_yyyyzz[k] = -g_z_0_yzzz_yyyyzz[k] * cd_y[k] + g_z_0_yzzz_yyyyyzz[k];

                g_z_0_yyzzz_yyyzzz[k] = -g_z_0_yzzz_yyyzzz[k] * cd_y[k] + g_z_0_yzzz_yyyyzzz[k];

                g_z_0_yyzzz_yyzzzz[k] = -g_z_0_yzzz_yyzzzz[k] * cd_y[k] + g_z_0_yzzz_yyyzzzz[k];

                g_z_0_yyzzz_yzzzzz[k] = -g_z_0_yzzz_yzzzzz[k] * cd_y[k] + g_z_0_yzzz_yyzzzzz[k];

                g_z_0_yyzzz_zzzzzz[k] = -g_z_0_yzzz_zzzzzz[k] * cd_y[k] + g_z_0_yzzz_yzzzzzz[k];
            }

            /// Set up 532-560 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 532);

            auto g_z_0_yzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 533);

            auto g_z_0_yzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 534);

            auto g_z_0_yzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 535);

            auto g_z_0_yzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 536);

            auto g_z_0_yzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 537);

            auto g_z_0_yzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 538);

            auto g_z_0_yzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 539);

            auto g_z_0_yzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 540);

            auto g_z_0_yzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 541);

            auto g_z_0_yzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 542);

            auto g_z_0_yzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 543);

            auto g_z_0_yzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 544);

            auto g_z_0_yzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 545);

            auto g_z_0_yzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 546);

            auto g_z_0_yzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 547);

            auto g_z_0_yzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 548);

            auto g_z_0_yzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 549);

            auto g_z_0_yzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 550);

            auto g_z_0_yzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 551);

            auto g_z_0_yzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 552);

            auto g_z_0_yzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 553);

            auto g_z_0_yzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 554);

            auto g_z_0_yzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 555);

            auto g_z_0_yzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 556);

            auto g_z_0_yzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 557);

            auto g_z_0_yzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 558);

            auto g_z_0_yzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 559);

            #pragma omp simd aligned(cd_y, g_z_0_yzzzz_xxxxxx, g_z_0_yzzzz_xxxxxy, g_z_0_yzzzz_xxxxxz, g_z_0_yzzzz_xxxxyy, g_z_0_yzzzz_xxxxyz, g_z_0_yzzzz_xxxxzz, g_z_0_yzzzz_xxxyyy, g_z_0_yzzzz_xxxyyz, g_z_0_yzzzz_xxxyzz, g_z_0_yzzzz_xxxzzz, g_z_0_yzzzz_xxyyyy, g_z_0_yzzzz_xxyyyz, g_z_0_yzzzz_xxyyzz, g_z_0_yzzzz_xxyzzz, g_z_0_yzzzz_xxzzzz, g_z_0_yzzzz_xyyyyy, g_z_0_yzzzz_xyyyyz, g_z_0_yzzzz_xyyyzz, g_z_0_yzzzz_xyyzzz, g_z_0_yzzzz_xyzzzz, g_z_0_yzzzz_xzzzzz, g_z_0_yzzzz_yyyyyy, g_z_0_yzzzz_yyyyyz, g_z_0_yzzzz_yyyyzz, g_z_0_yzzzz_yyyzzz, g_z_0_yzzzz_yyzzzz, g_z_0_yzzzz_yzzzzz, g_z_0_yzzzz_zzzzzz, g_z_0_zzzz_xxxxxx, g_z_0_zzzz_xxxxxxy, g_z_0_zzzz_xxxxxy, g_z_0_zzzz_xxxxxyy, g_z_0_zzzz_xxxxxyz, g_z_0_zzzz_xxxxxz, g_z_0_zzzz_xxxxyy, g_z_0_zzzz_xxxxyyy, g_z_0_zzzz_xxxxyyz, g_z_0_zzzz_xxxxyz, g_z_0_zzzz_xxxxyzz, g_z_0_zzzz_xxxxzz, g_z_0_zzzz_xxxyyy, g_z_0_zzzz_xxxyyyy, g_z_0_zzzz_xxxyyyz, g_z_0_zzzz_xxxyyz, g_z_0_zzzz_xxxyyzz, g_z_0_zzzz_xxxyzz, g_z_0_zzzz_xxxyzzz, g_z_0_zzzz_xxxzzz, g_z_0_zzzz_xxyyyy, g_z_0_zzzz_xxyyyyy, g_z_0_zzzz_xxyyyyz, g_z_0_zzzz_xxyyyz, g_z_0_zzzz_xxyyyzz, g_z_0_zzzz_xxyyzz, g_z_0_zzzz_xxyyzzz, g_z_0_zzzz_xxyzzz, g_z_0_zzzz_xxyzzzz, g_z_0_zzzz_xxzzzz, g_z_0_zzzz_xyyyyy, g_z_0_zzzz_xyyyyyy, g_z_0_zzzz_xyyyyyz, g_z_0_zzzz_xyyyyz, g_z_0_zzzz_xyyyyzz, g_z_0_zzzz_xyyyzz, g_z_0_zzzz_xyyyzzz, g_z_0_zzzz_xyyzzz, g_z_0_zzzz_xyyzzzz, g_z_0_zzzz_xyzzzz, g_z_0_zzzz_xyzzzzz, g_z_0_zzzz_xzzzzz, g_z_0_zzzz_yyyyyy, g_z_0_zzzz_yyyyyyy, g_z_0_zzzz_yyyyyyz, g_z_0_zzzz_yyyyyz, g_z_0_zzzz_yyyyyzz, g_z_0_zzzz_yyyyzz, g_z_0_zzzz_yyyyzzz, g_z_0_zzzz_yyyzzz, g_z_0_zzzz_yyyzzzz, g_z_0_zzzz_yyzzzz, g_z_0_zzzz_yyzzzzz, g_z_0_zzzz_yzzzzz, g_z_0_zzzz_yzzzzzz, g_z_0_zzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzz_xxxxxx[k] = -g_z_0_zzzz_xxxxxx[k] * cd_y[k] + g_z_0_zzzz_xxxxxxy[k];

                g_z_0_yzzzz_xxxxxy[k] = -g_z_0_zzzz_xxxxxy[k] * cd_y[k] + g_z_0_zzzz_xxxxxyy[k];

                g_z_0_yzzzz_xxxxxz[k] = -g_z_0_zzzz_xxxxxz[k] * cd_y[k] + g_z_0_zzzz_xxxxxyz[k];

                g_z_0_yzzzz_xxxxyy[k] = -g_z_0_zzzz_xxxxyy[k] * cd_y[k] + g_z_0_zzzz_xxxxyyy[k];

                g_z_0_yzzzz_xxxxyz[k] = -g_z_0_zzzz_xxxxyz[k] * cd_y[k] + g_z_0_zzzz_xxxxyyz[k];

                g_z_0_yzzzz_xxxxzz[k] = -g_z_0_zzzz_xxxxzz[k] * cd_y[k] + g_z_0_zzzz_xxxxyzz[k];

                g_z_0_yzzzz_xxxyyy[k] = -g_z_0_zzzz_xxxyyy[k] * cd_y[k] + g_z_0_zzzz_xxxyyyy[k];

                g_z_0_yzzzz_xxxyyz[k] = -g_z_0_zzzz_xxxyyz[k] * cd_y[k] + g_z_0_zzzz_xxxyyyz[k];

                g_z_0_yzzzz_xxxyzz[k] = -g_z_0_zzzz_xxxyzz[k] * cd_y[k] + g_z_0_zzzz_xxxyyzz[k];

                g_z_0_yzzzz_xxxzzz[k] = -g_z_0_zzzz_xxxzzz[k] * cd_y[k] + g_z_0_zzzz_xxxyzzz[k];

                g_z_0_yzzzz_xxyyyy[k] = -g_z_0_zzzz_xxyyyy[k] * cd_y[k] + g_z_0_zzzz_xxyyyyy[k];

                g_z_0_yzzzz_xxyyyz[k] = -g_z_0_zzzz_xxyyyz[k] * cd_y[k] + g_z_0_zzzz_xxyyyyz[k];

                g_z_0_yzzzz_xxyyzz[k] = -g_z_0_zzzz_xxyyzz[k] * cd_y[k] + g_z_0_zzzz_xxyyyzz[k];

                g_z_0_yzzzz_xxyzzz[k] = -g_z_0_zzzz_xxyzzz[k] * cd_y[k] + g_z_0_zzzz_xxyyzzz[k];

                g_z_0_yzzzz_xxzzzz[k] = -g_z_0_zzzz_xxzzzz[k] * cd_y[k] + g_z_0_zzzz_xxyzzzz[k];

                g_z_0_yzzzz_xyyyyy[k] = -g_z_0_zzzz_xyyyyy[k] * cd_y[k] + g_z_0_zzzz_xyyyyyy[k];

                g_z_0_yzzzz_xyyyyz[k] = -g_z_0_zzzz_xyyyyz[k] * cd_y[k] + g_z_0_zzzz_xyyyyyz[k];

                g_z_0_yzzzz_xyyyzz[k] = -g_z_0_zzzz_xyyyzz[k] * cd_y[k] + g_z_0_zzzz_xyyyyzz[k];

                g_z_0_yzzzz_xyyzzz[k] = -g_z_0_zzzz_xyyzzz[k] * cd_y[k] + g_z_0_zzzz_xyyyzzz[k];

                g_z_0_yzzzz_xyzzzz[k] = -g_z_0_zzzz_xyzzzz[k] * cd_y[k] + g_z_0_zzzz_xyyzzzz[k];

                g_z_0_yzzzz_xzzzzz[k] = -g_z_0_zzzz_xzzzzz[k] * cd_y[k] + g_z_0_zzzz_xyzzzzz[k];

                g_z_0_yzzzz_yyyyyy[k] = -g_z_0_zzzz_yyyyyy[k] * cd_y[k] + g_z_0_zzzz_yyyyyyy[k];

                g_z_0_yzzzz_yyyyyz[k] = -g_z_0_zzzz_yyyyyz[k] * cd_y[k] + g_z_0_zzzz_yyyyyyz[k];

                g_z_0_yzzzz_yyyyzz[k] = -g_z_0_zzzz_yyyyzz[k] * cd_y[k] + g_z_0_zzzz_yyyyyzz[k];

                g_z_0_yzzzz_yyyzzz[k] = -g_z_0_zzzz_yyyzzz[k] * cd_y[k] + g_z_0_zzzz_yyyyzzz[k];

                g_z_0_yzzzz_yyzzzz[k] = -g_z_0_zzzz_yyzzzz[k] * cd_y[k] + g_z_0_zzzz_yyyzzzz[k];

                g_z_0_yzzzz_yzzzzz[k] = -g_z_0_zzzz_yzzzzz[k] * cd_y[k] + g_z_0_zzzz_yyzzzzz[k];

                g_z_0_yzzzz_zzzzzz[k] = -g_z_0_zzzz_zzzzzz[k] * cd_y[k] + g_z_0_zzzz_yzzzzzz[k];
            }

            /// Set up 560-588 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 560);

            auto g_z_0_zzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 561);

            auto g_z_0_zzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 562);

            auto g_z_0_zzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 563);

            auto g_z_0_zzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 564);

            auto g_z_0_zzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 565);

            auto g_z_0_zzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 566);

            auto g_z_0_zzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 567);

            auto g_z_0_zzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 568);

            auto g_z_0_zzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 569);

            auto g_z_0_zzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 570);

            auto g_z_0_zzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 571);

            auto g_z_0_zzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 572);

            auto g_z_0_zzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 573);

            auto g_z_0_zzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 574);

            auto g_z_0_zzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 575);

            auto g_z_0_zzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 576);

            auto g_z_0_zzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 577);

            auto g_z_0_zzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 578);

            auto g_z_0_zzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 579);

            auto g_z_0_zzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 580);

            auto g_z_0_zzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 581);

            auto g_z_0_zzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 582);

            auto g_z_0_zzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 583);

            auto g_z_0_zzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 584);

            auto g_z_0_zzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 585);

            auto g_z_0_zzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 586);

            auto g_z_0_zzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 587);

            #pragma omp simd aligned(cd_z, g_z_0_zzzz_xxxxxx, g_z_0_zzzz_xxxxxxz, g_z_0_zzzz_xxxxxy, g_z_0_zzzz_xxxxxyz, g_z_0_zzzz_xxxxxz, g_z_0_zzzz_xxxxxzz, g_z_0_zzzz_xxxxyy, g_z_0_zzzz_xxxxyyz, g_z_0_zzzz_xxxxyz, g_z_0_zzzz_xxxxyzz, g_z_0_zzzz_xxxxzz, g_z_0_zzzz_xxxxzzz, g_z_0_zzzz_xxxyyy, g_z_0_zzzz_xxxyyyz, g_z_0_zzzz_xxxyyz, g_z_0_zzzz_xxxyyzz, g_z_0_zzzz_xxxyzz, g_z_0_zzzz_xxxyzzz, g_z_0_zzzz_xxxzzz, g_z_0_zzzz_xxxzzzz, g_z_0_zzzz_xxyyyy, g_z_0_zzzz_xxyyyyz, g_z_0_zzzz_xxyyyz, g_z_0_zzzz_xxyyyzz, g_z_0_zzzz_xxyyzz, g_z_0_zzzz_xxyyzzz, g_z_0_zzzz_xxyzzz, g_z_0_zzzz_xxyzzzz, g_z_0_zzzz_xxzzzz, g_z_0_zzzz_xxzzzzz, g_z_0_zzzz_xyyyyy, g_z_0_zzzz_xyyyyyz, g_z_0_zzzz_xyyyyz, g_z_0_zzzz_xyyyyzz, g_z_0_zzzz_xyyyzz, g_z_0_zzzz_xyyyzzz, g_z_0_zzzz_xyyzzz, g_z_0_zzzz_xyyzzzz, g_z_0_zzzz_xyzzzz, g_z_0_zzzz_xyzzzzz, g_z_0_zzzz_xzzzzz, g_z_0_zzzz_xzzzzzz, g_z_0_zzzz_yyyyyy, g_z_0_zzzz_yyyyyyz, g_z_0_zzzz_yyyyyz, g_z_0_zzzz_yyyyyzz, g_z_0_zzzz_yyyyzz, g_z_0_zzzz_yyyyzzz, g_z_0_zzzz_yyyzzz, g_z_0_zzzz_yyyzzzz, g_z_0_zzzz_yyzzzz, g_z_0_zzzz_yyzzzzz, g_z_0_zzzz_yzzzzz, g_z_0_zzzz_yzzzzzz, g_z_0_zzzz_zzzzzz, g_z_0_zzzz_zzzzzzz, g_z_0_zzzzz_xxxxxx, g_z_0_zzzzz_xxxxxy, g_z_0_zzzzz_xxxxxz, g_z_0_zzzzz_xxxxyy, g_z_0_zzzzz_xxxxyz, g_z_0_zzzzz_xxxxzz, g_z_0_zzzzz_xxxyyy, g_z_0_zzzzz_xxxyyz, g_z_0_zzzzz_xxxyzz, g_z_0_zzzzz_xxxzzz, g_z_0_zzzzz_xxyyyy, g_z_0_zzzzz_xxyyyz, g_z_0_zzzzz_xxyyzz, g_z_0_zzzzz_xxyzzz, g_z_0_zzzzz_xxzzzz, g_z_0_zzzzz_xyyyyy, g_z_0_zzzzz_xyyyyz, g_z_0_zzzzz_xyyyzz, g_z_0_zzzzz_xyyzzz, g_z_0_zzzzz_xyzzzz, g_z_0_zzzzz_xzzzzz, g_z_0_zzzzz_yyyyyy, g_z_0_zzzzz_yyyyyz, g_z_0_zzzzz_yyyyzz, g_z_0_zzzzz_yyyzzz, g_z_0_zzzzz_yyzzzz, g_z_0_zzzzz_yzzzzz, g_z_0_zzzzz_zzzzzz, g_zzzz_xxxxxx, g_zzzz_xxxxxy, g_zzzz_xxxxxz, g_zzzz_xxxxyy, g_zzzz_xxxxyz, g_zzzz_xxxxzz, g_zzzz_xxxyyy, g_zzzz_xxxyyz, g_zzzz_xxxyzz, g_zzzz_xxxzzz, g_zzzz_xxyyyy, g_zzzz_xxyyyz, g_zzzz_xxyyzz, g_zzzz_xxyzzz, g_zzzz_xxzzzz, g_zzzz_xyyyyy, g_zzzz_xyyyyz, g_zzzz_xyyyzz, g_zzzz_xyyzzz, g_zzzz_xyzzzz, g_zzzz_xzzzzz, g_zzzz_yyyyyy, g_zzzz_yyyyyz, g_zzzz_yyyyzz, g_zzzz_yyyzzz, g_zzzz_yyzzzz, g_zzzz_yzzzzz, g_zzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzz_xxxxxx[k] = -g_zzzz_xxxxxx[k] - g_z_0_zzzz_xxxxxx[k] * cd_z[k] + g_z_0_zzzz_xxxxxxz[k];

                g_z_0_zzzzz_xxxxxy[k] = -g_zzzz_xxxxxy[k] - g_z_0_zzzz_xxxxxy[k] * cd_z[k] + g_z_0_zzzz_xxxxxyz[k];

                g_z_0_zzzzz_xxxxxz[k] = -g_zzzz_xxxxxz[k] - g_z_0_zzzz_xxxxxz[k] * cd_z[k] + g_z_0_zzzz_xxxxxzz[k];

                g_z_0_zzzzz_xxxxyy[k] = -g_zzzz_xxxxyy[k] - g_z_0_zzzz_xxxxyy[k] * cd_z[k] + g_z_0_zzzz_xxxxyyz[k];

                g_z_0_zzzzz_xxxxyz[k] = -g_zzzz_xxxxyz[k] - g_z_0_zzzz_xxxxyz[k] * cd_z[k] + g_z_0_zzzz_xxxxyzz[k];

                g_z_0_zzzzz_xxxxzz[k] = -g_zzzz_xxxxzz[k] - g_z_0_zzzz_xxxxzz[k] * cd_z[k] + g_z_0_zzzz_xxxxzzz[k];

                g_z_0_zzzzz_xxxyyy[k] = -g_zzzz_xxxyyy[k] - g_z_0_zzzz_xxxyyy[k] * cd_z[k] + g_z_0_zzzz_xxxyyyz[k];

                g_z_0_zzzzz_xxxyyz[k] = -g_zzzz_xxxyyz[k] - g_z_0_zzzz_xxxyyz[k] * cd_z[k] + g_z_0_zzzz_xxxyyzz[k];

                g_z_0_zzzzz_xxxyzz[k] = -g_zzzz_xxxyzz[k] - g_z_0_zzzz_xxxyzz[k] * cd_z[k] + g_z_0_zzzz_xxxyzzz[k];

                g_z_0_zzzzz_xxxzzz[k] = -g_zzzz_xxxzzz[k] - g_z_0_zzzz_xxxzzz[k] * cd_z[k] + g_z_0_zzzz_xxxzzzz[k];

                g_z_0_zzzzz_xxyyyy[k] = -g_zzzz_xxyyyy[k] - g_z_0_zzzz_xxyyyy[k] * cd_z[k] + g_z_0_zzzz_xxyyyyz[k];

                g_z_0_zzzzz_xxyyyz[k] = -g_zzzz_xxyyyz[k] - g_z_0_zzzz_xxyyyz[k] * cd_z[k] + g_z_0_zzzz_xxyyyzz[k];

                g_z_0_zzzzz_xxyyzz[k] = -g_zzzz_xxyyzz[k] - g_z_0_zzzz_xxyyzz[k] * cd_z[k] + g_z_0_zzzz_xxyyzzz[k];

                g_z_0_zzzzz_xxyzzz[k] = -g_zzzz_xxyzzz[k] - g_z_0_zzzz_xxyzzz[k] * cd_z[k] + g_z_0_zzzz_xxyzzzz[k];

                g_z_0_zzzzz_xxzzzz[k] = -g_zzzz_xxzzzz[k] - g_z_0_zzzz_xxzzzz[k] * cd_z[k] + g_z_0_zzzz_xxzzzzz[k];

                g_z_0_zzzzz_xyyyyy[k] = -g_zzzz_xyyyyy[k] - g_z_0_zzzz_xyyyyy[k] * cd_z[k] + g_z_0_zzzz_xyyyyyz[k];

                g_z_0_zzzzz_xyyyyz[k] = -g_zzzz_xyyyyz[k] - g_z_0_zzzz_xyyyyz[k] * cd_z[k] + g_z_0_zzzz_xyyyyzz[k];

                g_z_0_zzzzz_xyyyzz[k] = -g_zzzz_xyyyzz[k] - g_z_0_zzzz_xyyyzz[k] * cd_z[k] + g_z_0_zzzz_xyyyzzz[k];

                g_z_0_zzzzz_xyyzzz[k] = -g_zzzz_xyyzzz[k] - g_z_0_zzzz_xyyzzz[k] * cd_z[k] + g_z_0_zzzz_xyyzzzz[k];

                g_z_0_zzzzz_xyzzzz[k] = -g_zzzz_xyzzzz[k] - g_z_0_zzzz_xyzzzz[k] * cd_z[k] + g_z_0_zzzz_xyzzzzz[k];

                g_z_0_zzzzz_xzzzzz[k] = -g_zzzz_xzzzzz[k] - g_z_0_zzzz_xzzzzz[k] * cd_z[k] + g_z_0_zzzz_xzzzzzz[k];

                g_z_0_zzzzz_yyyyyy[k] = -g_zzzz_yyyyyy[k] - g_z_0_zzzz_yyyyyy[k] * cd_z[k] + g_z_0_zzzz_yyyyyyz[k];

                g_z_0_zzzzz_yyyyyz[k] = -g_zzzz_yyyyyz[k] - g_z_0_zzzz_yyyyyz[k] * cd_z[k] + g_z_0_zzzz_yyyyyzz[k];

                g_z_0_zzzzz_yyyyzz[k] = -g_zzzz_yyyyzz[k] - g_z_0_zzzz_yyyyzz[k] * cd_z[k] + g_z_0_zzzz_yyyyzzz[k];

                g_z_0_zzzzz_yyyzzz[k] = -g_zzzz_yyyzzz[k] - g_z_0_zzzz_yyyzzz[k] * cd_z[k] + g_z_0_zzzz_yyyzzzz[k];

                g_z_0_zzzzz_yyzzzz[k] = -g_zzzz_yyzzzz[k] - g_z_0_zzzz_yyzzzz[k] * cd_z[k] + g_z_0_zzzz_yyzzzzz[k];

                g_z_0_zzzzz_yzzzzz[k] = -g_zzzz_yzzzzz[k] - g_z_0_zzzz_yzzzzz[k] * cd_z[k] + g_z_0_zzzz_yzzzzzz[k];

                g_z_0_zzzzz_zzzzzz[k] = -g_zzzz_zzzzzz[k] - g_z_0_zzzz_zzzzzz[k] * cd_z[k] + g_z_0_zzzz_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

