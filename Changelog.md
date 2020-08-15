# Changelog

## Half-Life: Static Friction

### 2020-08-15
* Removed `./GameSource/hlsf/sound/music/Death.flac`
  * Was copyrighted.

### 2012-01-28

* Added slow-motion capability
  * Accessed through func_slowmo

### 2012-01-26

* Fixed Engine Access Violation ocurring in hlsf.dll
  * Visual Studio's optimization was doing something
    wacky to the code, disable it until I can figure
    out why.

### 2012-1-23

* Removed item_gravsuit
  * Replaced (see below)

* Added new input type system
  * New Entity: input_settype

### 2012-1-22

* Added HEV Barney
  * Set by a keyvalue in monster_barney
* Screen no longer stays red after death in multiplayer

* Added item_gravsuit
  * Does some funky stuff to your view controls.

### 2012-1-21

* Added env_deadplayer
  * This is created when the player dies
* Reenabled multiplayer support.
  * vsys_* entities don't function in multiplayer
  * Screen stays faded after respawn in multiplayer
* Added monster_kate and monster_barniel and their respective "dead" entities
  * Uses Barney's AI
* Added func_detail entity.
  * Alias for func_wall
* Added everyone's old friend: Otis!
  * monster_otis
* Friendly NPC's now try to talk to the player.
* Added monster_gina
* Added monster_colette

### 2012-1-20

* Added dead model(s) to death sequence
  * Selects appropriate model depending on whether or not
    the player has item_supersuit
  * Hidden actual player model on death ("models/player.mdl")
* Headcrabs now spawn at the head of the zombie, not the (chest?)
* Removed vsys_release.
  * Pointless.
* Temporarily disabled HEV Introduction.
  * Causes game to crash.

### 2012-1-19

* Added Death fade-to-red
  * Displays 'Mission Failed'
  * Changes view to thirdperson
  
* All commands that attempt to Play a CD Audio track are redirected to VSys

* Func_recharge now only allows recharge with item_supersuit

* Weapon_handgrenade only gives 1 by default instead of 5

* Renamed vsys_setup to vsys_object
  * Now toggles sound on/off
  * Other vsys_* entities are still functional
  
* Headcrabs now jump off zombies
  * If headshot, headcrab doesn't spawn

* Updated Zombie Model

### 2012-1-18 - Beta 0.01

* Initial setup with Spirit of Half-Life 1.8

* Added VSys Audio system.
  * Supports up to 32 channels
  
* Added VSys Management Entities
  * vsys_changevol
  * vsys_toggle
  * vsys_pause
  * vsys_resume
  * vsys_stop
  * vsys_play
  * vsys_release
  * vsys_setup

* item_suit no longer allows use of the flashlight, armor and batteries

* Added item_supersuit
  * Allows use of flashlight, armor and batteries
  
* Added Death music
  
* New monster_cleansuit_scientist
  * Same as monster_scientist except with a cleansuit
  
- Disabled multiplayer support
  * Crashed on load

## Half-Life: Sci-Fi

Precursor to Static Friction. Essentially a testing ground.

### 2008-02-26

* Added random ammo amounts to all weapons with clips, picking up a weapon (not ammo) gives a random number
* Can also be set via editor, NOTE that when setting in editor, must be within normal clip size or it doesn't work (over-ruled by code)
* New soldier zombie, select via variable in monster_zombie (automatically has 2x normal health)
* New prop_hevcharger, new interactive model HEV charger
* New prop_healthcharger, new interactive model health charger
* New monster_baby_gargantua, similar to normal garg only not immune to bullets
* New monster_male_assassin, uses 9mmAR and sniper rifle
* New monster_pitdrone
* Not fully functional, does not shoot yet
* New monster_shocktrooper
* Working but code needs tweaking - doesn't throw grenades
* New LAW (Light Anti-Tank Weapon) weapon, disposable, single shot rocket launcher, best on human targets
* New monster_human_grunt_law, Grunt with law weapon, fires once then switches to pistol
* Replaced RPG with new model
* New .aur files (see new catelouge.txt in aurora folder for reference)
* Improved spark sound effects and variety
* Modified player crouch behavior
* Modified tripmine beam
* Tweaked grenade behavior and physics, more realistic
* env_spark emits light by default, has an option to not
* monster_scientist_dead can now use cleansuit model, choose via a variable in editor
* All music now placed in music folder in root of mod, fmod responds better
* Changed name of mod to "Half-Life: Static Friction"
* Zombie head changed to more realistic version
* M40 Sniper Rifle now zooms twice and resets properly
* Fixed error in keyvalue that set wether Barney uses HEV suit or not
* Fixed bouncy crouch bug (thanks to Cale 'Mazor' Dunlap)
* Fixed bug where M40's zoomed crosshair did not disappear when reloading when zoomed in
* Barney's make proper shooting sound
* titles.txt file was missing vital entries, replaced with Half-Life original


### 2008-2-23 - Beta 0.04
* Replaced old decals with new better-looking ones
* Replaced 9mm Handgun with new Glock18 model and sounds
* Replaced 9mm Sub-Machinegun with new model and sounds
* Replaced Crossbow with new model and sounds
* Replaced Shotgun with new model and sounds
* Replaced Snark with new model and sounds
* Replaced 357 Magnum with new model and sounds
* Replaced Hornet Gun with new model and sounds
* Replaced Egon with new model
* Replaced Gauss with new model
* Replaced Grenade with new model and sounds
* New Gina character from various Half-Life Episodes, uses 357 Magnum
* Added recoil effect to all applicable weapons
* New M40 Rifle
* New Particle (.aur) effect, "burning_smoke.aur". Use to create smoke effects in maps
* Included more Hi-Def models
* New Zombie blood effects, shooting body spawns red blood, shooting headcrab spawns yellow blood
* New debris effects added to explosions
* Improved explosion effects
* Added slot 6
* Modified crowbar model (again!)
* 9mm Handgun is now properly semi-automatic, holding fire doesn't keep shooting
* This means, however, that you can fire as fast as you can pull the trigger making secondary fire pointless
* 9mm Handgun does less damage to reflect new firing system
* Tweaked floating code
* Fixed Barney zombie model, headcrab was sitting too high
* Hud colour adjusted
* It is now a rather appealing glowing cyan, I was aiming for white but got this instead!
* Monsters no longer attack dead enemies
* Monsters on ground are no longer solid
* Fixed Zombie hitbox code, shooting different body parts now correctly adjusts inflicted damage
* Fixed a conflict between floating code and zombie spawn crab feature
* Fixed spawn position of Zombie Headcrabs to head instead of chest
* Headcrabs now "jump" off zombies instead of falling to the floor

### 2008-2-9 - Beta 0.03

* Barnacle now plays idle sounds and a variety of pulling, chewing and digesting sounds
* Alternative version of Barney, wears a HEV Suit, can choose via a variable in the monster_barney
* New trigger_timeadjust, adjusts game speed to allow for super fast or slowmo gameplay
* Implemented Cale 'Mazor' Dunlap's view roll when strafing code
* New dead body floating code implemented in all monsters
* Code is stable but there are bugs that need fixing
* Bullet casings now lie on the ground for 60 seconds instead of 2.5 seconds
* Colette character finished, speech converted to HL1 format, uses 357 Magnum
* Animation error on Fast Headcrab model fixed


### 2008-1-30 - Beta 0.02

* New Fast Headcrab monster, operates same as headcrab
* Added "autoexec.cfg" to display text when console appears
* Included Hi-Definition monsters
* Alternative version of scientist, wears a cleansuit, can choose via a variable in the monster_scientist
* Added Zombie death sounds
* Headcrabs now jump off dead zombies
* Alternative version of zombie, barney zombie, can choose via a variable in the monster_zombie
* New Colette Green character from Half-Life Decay
* Colette is not finished so does not appear in tester map
* Changed old zombie sounds into new HL2 zombie sounds
* Modified crowbar sounds and models
* Changed mod dll name so it says "Dll loaded for mod: Half-Life: SciFi" instead of "Spirit of Half-Life"
* Gordon's muzzleflash moved to correct position


### 2008-1-23 - Beta 0.01a

* Added Blue-Shift Rain Sprites
* Added Custom Scientist Voices
* Added Rosenberg.wad
* Added 2 Songs of Soundtrack
* Added Custom Scientist Voices in sentences.txt
* Fixed liblist.gam to show "Half-Life: Sci-Fi" instead of "Spirit of Half-Life"


### 2008-??-??

* Initial setup with Spirit of Half-Life 1.5alpha4
* Included origial menu Splash Screens and buttons
* New Gordon Freeman character, uses shotgun, doesn't speak
* New Barniel character, uses new Scorpian Uzi, new speech dialog